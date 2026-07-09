#!/usr/bin/env python3

"""
Parse KofamScan ribosomal annotations and convert protein-level hits into
gene-level ribosomal reference annotations.

KofamScan searches amino-acid sequences and therefore reports annotations
using protein identifiers. Downstream codon-usage analyses operate on genes
and coding DNA sequences. This script connects these two identifier spaces
using a gene-to-protein mapping table.

Only KofamScan hits marked as significant are retained. In the
``detail-tsv`` output, a significant hit is indicated by a leading asterisk
and represents a score passing the KO-specific KOfam threshold.

The script verifies that every reported KO belongs to the restricted
ribosomal KO reference and that every significant protein can be mapped to
a gene identifier.

Multiple significant KO assignments are preserved rather than collapsed
arbitrarily. This allows ambiguous annotations and paralogous ribosomal
proteins to be identified during quality control.

The script produces:

``ribosome_significant_hits.tsv``
    One row per significant protein-KO hit, including protein ID, gene ID,
    KO identifier, score, threshold, E-value and functional definition.

``ribosome_gene_annotations.tsv``
    Gene-level aggregation of the significant KofamScan hits. The table
    contains all assigned KO identifiers and flags genes with multiple KO
    assignments.

``ribosomal_reference_gene_ids.txt``
    A plain-text list of gene identifiers classified as cytosolic
    ribosomal genes.

``ribosome_qc.tsv``
    Per-genome summary containing the number of proteins, detected
    ribosomal genes, observed KO identifiers, KO coverage, ambiguous hits
    and unmapped proteins.

Snakemake inputs
----------------
detail
    KofamScan ``detail-tsv`` output generated using the restricted
    ribosomal HAL profile set.
proteins
    Protein FASTA file used as the KofamScan query.
mapping
    Tab-separated mapping between protein identifiers and gene identifiers.
ko_table
    Table of allowed cytosolic eukaryotic ribosomal KO identifiers.

Snakemake outputs
-----------------
hits
    Significant protein-level KofamScan annotations.
genes
    Aggregated gene-level ribosomal annotations.
gene_ids
    Gene identifiers selected for the ribosomal reference set.
qc
    Per-sample annotation quality-control summary.

Snakemake parameters
--------------------
gene_column
    Name of the gene identifier column in the mapping table.
protein_column
    Name of the protein identifier column in the mapping table.

Why this script is required
---------------------------
Raw KofamScan output cannot be used directly to select coding sequences.
It contains protein identifiers, may contain non-significant hits and may
assign multiple KO identifiers to one sequence. This script applies the
Kofam thresholds, validates the restricted KO scope, maps proteins to genes
and creates an auditable list of ribosomal reference genes.
"""

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Any, TextIO

from Bio import SeqIO

KO_PATTERN = re.compile(r"^K\d{5}$")


def open_text(path: Path) -> TextIO:
    """Open plain-text or gzip-compressed input."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_reference_kos(path: Path) -> set[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "ko" not in reader.fieldnames:
            raise ValueError(f"Reference KO table must contain a 'ko' column: {path}")
        kos = {row["ko"].strip() for row in reader if row.get("ko", "").strip()}

    invalid = sorted(ko for ko in kos if not KO_PATTERN.fullmatch(ko))
    if invalid:
        raise ValueError(f"Invalid KO identifiers in {path}: {', '.join(invalid[:20])}")
    if not kos:
        raise ValueError(f"No KO identifiers found in {path}")
    return kos


def read_gene_protein_map(
    path: Path,
    gene_column: str,
    protein_column: str,
    cds_column: str,
) -> dict[str, dict[str, str]]:
    """Return protein_id -> {gene_id, cds_id} with one-to-one validation."""
    mapping: dict[str, dict[str, str]] = {}

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")

        required = {gene_column, protein_column, cds_column}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Mapping table {path} lacks required columns: "
                f"{', '.join(sorted(missing))}"
            )

        for row in reader:
            gene_id = row.get(gene_column, "").strip()
            protein_id = row.get(protein_column, "").strip()
            cds_id = row.get(cds_column, "").strip()

            if not gene_id or not protein_id or not cds_id:
                continue

            current = {
                "gene_id": gene_id,
                "cds_id": cds_id,
            }

            previous = mapping.get(protein_id)
            if previous is not None and previous != current:
                raise ValueError(
                    f"Protein {protein_id} maps to multiple gene/CDS records: "
                    f"{previous} and {current}"
                )

            mapping[protein_id] = current

    if not mapping:
        raise ValueError(f"No complete protein-to-gene/CDS mappings found in {path}")

    return mapping


def filter_mapping_to_proteome(
    protein_to_ids: dict[str, dict[str, str]],
    protein_ids: set[str],
    max_missing_count: int,
    max_missing_fraction: float,
) -> tuple[dict[str, dict[str, str]], int, float]:
    """Remove mappings absent from the query proteome unless mismatch is severe."""
    mapped_protein_ids = set(protein_to_ids)
    missing_from_proteome = sorted(mapped_protein_ids - protein_ids)

    missing_count = len(missing_from_proteome)
    missing_fraction = missing_count / max(len(mapped_protein_ids), 1)

    if missing_count > 0:
        preview = ", ".join(missing_from_proteome[:20])
        message = (
            f"{missing_count} mapped protein IDs are absent from the query proteome. "
            f"Missing fraction: {missing_fraction:.6f}. "
            f"Examples: {preview}"
        )

        if (
            missing_count > max_missing_count
            and missing_fraction > max_missing_fraction
        ):
            raise RuntimeError(
                message
                + (
                    f" Exceeds allowed thresholds: "
                    f"max_missing_proteome_count={max_missing_count}, "
                    f"max_missing_proteome_fraction={max_missing_fraction}."
                )
            )

        print("WARNING:", message)

    filtered_mapping = {
        protein_id: ids
        for protein_id, ids in protein_to_ids.items()
        if protein_id in protein_ids
    }

    if not filtered_mapping:
        raise RuntimeError(
            "No protein-to-gene/CDS mappings remained after filtering to "
            "protein IDs present in the query proteome."
        )

    return filtered_mapping, missing_count, missing_fraction


def read_unique_fasta_ids(path: Path) -> set[str]:
    ids: set[str] = set()
    duplicates: set[str] = set()
    with open_text(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in ids:
                duplicates.add(record.id)
            ids.add(record.id)

    if duplicates:
        raise ValueError(
            f"Duplicate protein FASTA identifiers in {path}: "
            + ", ".join(sorted(duplicates)[:20])
        )
    if not ids:
        raise ValueError(f"No protein records found in {path}")
    return ids


def parse_hit(line: str) -> dict[str, Any] | None:
    """Parse one KofamScan detail-tsv row."""
    stripped = line.rstrip("\n")
    if not stripped.strip() or stripped.lstrip().startswith("#"):
        return None

    stripped = stripped.lstrip()
    significant = stripped.startswith("*")
    if significant:
        stripped = stripped[1:].lstrip()

    fields = stripped.split("\t", maxsplit=5)
    if len(fields) != 6:
        fields = re.split(r"\s+", stripped, maxsplit=5)
    if len(fields) != 6:
        raise ValueError(f"Cannot parse KofamScan line: {line!r}")

    protein_id, ko_id, threshold, score, evalue, definition = fields
    if not KO_PATTERN.fullmatch(ko_id):
        raise ValueError(f"Invalid KO identifier in KofamScan output: {ko_id}")

    try:
        score_value = float(score)
        evalue_value = float(evalue)
        threshold_value = float(threshold)
    except ValueError as exc:
        raise ValueError(f"Invalid numeric value in KofamScan line: {line!r}") from exc

    return {
        "significant": significant,
        "protein_id": protein_id,
        "ko": ko_id,
        "threshold": threshold_value,
        "score": score_value,
        "evalue": evalue_value,
        "definition": definition.strip().strip('"'),
    }


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def write_id_list(path: Path, values: set[str] | list[str]) -> None:
    values_sorted = sorted(set(values))
    if not values_sorted:
        raise RuntimeError(f"Refusing to write an empty identifier list: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        "".join(f"{value}\n" for value in values_sorted), encoding="utf-8"
    )
    temporary.replace(path)


def main() -> None:
    detail_path = Path(snakemake.input.detail)
    mapping_path = Path(snakemake.input.mapping)
    protein_fasta = Path(snakemake.input.proteins)
    ko_table = Path(snakemake.input.ko_table)

    hits_output = Path(snakemake.output.hits)
    genes_output = Path(snakemake.output.genes)
    gene_ids_output = Path(snakemake.output.gene_ids)
    cds_ids_output = Path(snakemake.output.cds_ids)
    qc_output = Path(snakemake.output.qc)

    sample = str(snakemake.wildcards.sample)
    gene_column = str(snakemake.params.gene_column)
    protein_column = str(snakemake.params.protein_column)
    cds_column = str(snakemake.params.cds_column)

    allowed_kos = read_reference_kos(ko_table)

    protein_to_ids = read_gene_protein_map(
        mapping_path,
        gene_column=gene_column,
        protein_column=protein_column,
        cds_column=cds_column,
    )

    protein_ids = read_unique_fasta_ids(protein_fasta)

    max_missing_count = int(
        getattr(
            snakemake.params,
            "max_missing_proteome_count",
            200,
        )
    )

    max_missing_fraction = float(
        getattr(
            snakemake.params,
            "max_missing_proteome_fraction",
            0.05,
        )
    )

    n_mapping_missing_from_proteome_total = len(set(protein_to_ids) - protein_ids)
    n_mapping_proteins_before_filtering = len(protein_to_ids)

    protein_to_ids, n_mapping_missing_from_proteome, mapping_missing_fraction = (
        filter_mapping_to_proteome(
            protein_to_ids=protein_to_ids,
            protein_ids=protein_ids,
            max_missing_count=max_missing_count,
            max_missing_fraction=max_missing_fraction,
        )
    )

    n_mapping_proteins_after_filtering = len(protein_to_ids)

    significant_hits: list[dict[str, Any]] = []
    unmapped_proteins: set[str] = set()
    unexpected_proteins: set[str] = set()

    with detail_path.open(encoding="utf-8") as handle:
        for line in handle:
            hit = parse_hit(line)
            if hit is None or not hit["significant"]:
                continue
            if hit["ko"] not in allowed_kos:
                raise RuntimeError(
                    f"Unexpected KO {hit['ko']} in restricted ribosomal output"
                )
            if hit["protein_id"] not in protein_ids:
                unexpected_proteins.add(hit["protein_id"])
                continue

            mapped = protein_to_ids.get(hit["protein_id"])
            if mapped is None:
                unmapped_proteins.add(hit["protein_id"])
                continue

            significant_hits.append(
                {
                    "sample": sample,
                    "gene_id": mapped["gene_id"],
                    "cds_id": mapped["cds_id"],
                    **hit,
                }
            )

    if unexpected_proteins:
        preview = ", ".join(sorted(unexpected_proteins)[:20])
        raise RuntimeError(
            f"KofamScan reported {len(unexpected_proteins)} protein IDs absent "
            f"from the query FASTA. Examples: {preview}"
        )
    if unmapped_proteins:
        preview = ", ".join(sorted(unmapped_proteins)[:20])
        raise RuntimeError(
            f"{len(unmapped_proteins)} significant proteins could not be mapped "
            f"to gene/CDS IDs. Examples: {preview}"
        )
    if not significant_hits:
        raise RuntimeError(
            f"No significant ribosomal KofamScan hits were found for sample {sample}"
        )

    # Remove exact duplicate rows while retaining genuine multiple KO assignments.
    unique_hits = {
        (
            row["gene_id"],
            row["cds_id"],
            row["protein_id"],
            row["ko"],
            row["threshold"],
            row["score"],
            row["evalue"],
            row["definition"],
        ): row
        for row in significant_hits
    }
    significant_hits = sorted(
        unique_hits.values(),
        key=lambda row: (row["gene_id"], row["protein_id"], row["ko"]),
    )

    write_tsv(
        hits_output,
        significant_hits,
        [
            "sample",
            "gene_id",
            "cds_id",
            "protein_id",
            "ko",
            "threshold",
            "score",
            "evalue",
            "definition",
        ],
    )

    hits_by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    hits_by_protein: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for hit in significant_hits:
        hits_by_gene[hit["gene_id"]].append(hit)
        hits_by_protein[hit["protein_id"]].append(hit)

    gene_rows: list[dict[str, Any]] = []
    for gene_id, hits in sorted(hits_by_gene.items()):
        ko_ids = sorted({hit["ko"] for hit in hits})
        proteins = sorted({hit["protein_id"] for hit in hits})
        cds_ids = sorted({hit["cds_id"] for hit in hits})
        gene_rows.append(
            {
                "sample": sample,
                "gene_id": gene_id,
                "protein_ids": ";".join(proteins),
                "cds_ids": ";".join(cds_ids),
                "kegg_ko": ";".join(ko_ids),
                "n_significant_hits": len(hits),
                "n_unique_kos": len(ko_ids),
                "ambiguous_ko_assignment": len(ko_ids) > 1,
                "ribosomal_reference": True,
            }
        )

    write_tsv(
        genes_output,
        gene_rows,
        [
            "sample",
            "gene_id",
            "protein_ids",
            "cds_ids",
            "kegg_ko",
            "n_significant_hits",
            "n_unique_kos",
            "ambiguous_ko_assignment",
            "ribosomal_reference",
        ],
    )

    reference_gene_ids = set(hits_by_gene)
    reference_cds_ids = {hit["cds_id"] for hit in significant_hits}
    write_id_list(gene_ids_output, reference_gene_ids)
    write_id_list(cds_ids_output, reference_cds_ids)

    observed_kos = {hit["ko"] for hit in significant_hits}
    multi_ko_proteins = sum(
        len({hit["ko"] for hit in hits}) > 1 for hits in hits_by_protein.values()
    )
    qc_rows = [
        {
            "sample": sample,
            "n_proteins_total": len(protein_ids),
            "n_reference_kos": len(allowed_kos),
            "n_observed_kos": len(observed_kos),
            "ko_coverage": len(observed_kos) / len(allowed_kos),
            "n_significant_hits": len(significant_hits),
            "n_ribosomal_proteins": len(hits_by_protein),
            "n_ribosomal_genes": len(reference_gene_ids),
            "n_reference_cds": len(reference_cds_ids),
            "n_multi_ko_proteins": multi_ko_proteins,
            "n_mapping_proteins_before_filtering": n_mapping_proteins_before_filtering,
            "n_mapping_proteins_after_filtering": n_mapping_proteins_after_filtering,
            "n_mapping_missing_from_proteome": n_mapping_missing_from_proteome,
            "mapping_missing_from_proteome_fraction": mapping_missing_fraction,
            "n_unmapped_significant_proteins": 0,
        }
    ]
    write_tsv(
        qc_output,
        qc_rows,
        [
            "sample",
            "n_proteins_total",
            "n_mapping_proteins_before_filtering",
            "n_mapping_proteins_after_filtering",
            "n_mapping_missing_from_proteome",
            "mapping_missing_from_proteome_fraction",
            "n_reference_kos",
            "n_observed_kos",
            "ko_coverage",
            "n_significant_hits",
            "n_ribosomal_proteins",
            "n_ribosomal_genes",
            "n_reference_cds",
            "n_multi_ko_proteins",
            "n_unmapped_significant_proteins",
        ],
    )


if __name__ == "__main__":
    main()
