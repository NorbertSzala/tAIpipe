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
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

from Bio import SeqIO

KO_PATTERN = re.compile(r"^K\d{5}$")


def read_reference_kos(path: Path) -> set[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["ko"].strip() for row in reader}


def read_gene_protein_map(
    path: Path,
    gene_column: str,
    protein_column: str,
) -> dict[str, str]:
    mapping: dict[str, str] = {}

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}")

        required = {gene_column, protein_column}
        missing = required - set(reader.fieldnames)

        if missing:
            raise ValueError(
                f"Missing columns in {path}: " + ", ".join(sorted(missing))
            )

        for row in reader:
            gene_id = row[gene_column].strip()
            protein_id = row[protein_column].strip()

            if not gene_id or not protein_id:
                continue

            previous = mapping.get(protein_id)

            if previous is not None and previous != gene_id:
                raise ValueError(
                    f"Protein {protein_id} maps to multiple genes: "
                    f"{previous}, {gene_id}"
                )

            mapping[protein_id] = gene_id

    return mapping


def parse_hit(line: str) -> dict[str, Any] | None:
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
        raise ValueError(f"Invalid KO identifier: {ko_id}")

    return {
        "significant": significant,
        "protein_id": protein_id,
        "ko": ko_id,
        "threshold": threshold,
        "score": float(score),
        "evalue": float(evalue),
        "definition": definition.strip().strip('"'),
    }


def write_tsv(
    path: Path,
    rows: list[dict[str, Any]],
    fieldnames: list[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    detail_path = Path(snakemake.input.detail)
    mapping_path = Path(snakemake.input.mapping)
    protein_fasta = Path(snakemake.input.proteins)
    ko_table = Path(snakemake.input.ko_table)

    hits_output = Path(snakemake.output.hits)
    genes_output = Path(snakemake.output.genes)
    ids_output = Path(snakemake.output.gene_ids)
    qc_output = Path(snakemake.output.qc)

    sample = str(snakemake.wildcards.sample)
    gene_column = str(snakemake.params.gene_column)
    protein_column = str(snakemake.params.protein_column)

    allowed_kos = read_reference_kos(ko_table)

    protein_to_gene = read_gene_protein_map(
        mapping_path,
        gene_column=gene_column,
        protein_column=protein_column,
    )

    protein_ids = {record.id for record in SeqIO.parse(protein_fasta, "fasta")}

    significant_hits: list[dict[str, Any]] = []
    unmapped_proteins: set[str] = set()

    with detail_path.open(encoding="utf-8") as handle:
        for line in handle:
            hit = parse_hit(line)

            if hit is None or not hit["significant"]:
                continue

            if hit["ko"] not in allowed_kos:
                raise RuntimeError(
                    f"Unexpected KO {hit['ko']} in restricted "
                    "ribosomal KofamScan output"
                )

            gene_id = protein_to_gene.get(hit["protein_id"])

            if gene_id is None:
                unmapped_proteins.add(hit["protein_id"])
                continue

            significant_hits.append(
                {
                    "sample": sample,
                    "gene_id": gene_id,
                    **hit,
                }
            )

    if unmapped_proteins:
        examples = ", ".join(sorted(unmapped_proteins)[:20])

        raise RuntimeError(
            f"{len(unmapped_proteins)} significant proteins could "
            f"not be mapped to genes. Examples: {examples}"
        )

    significant_hits.sort(
        key=lambda row: (
            row["gene_id"],
            row["protein_id"],
            row["ko"],
        )
    )

    write_tsv(
        hits_output,
        significant_hits,
        [
            "sample",
            "gene_id",
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
        protein_hits = sorted({hit["protein_id"] for hit in hits})

        gene_rows.append(
            {
                "sample": sample,
                "gene_id": gene_id,
                "protein_ids": ";".join(protein_hits),
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
            "kegg_ko",
            "n_significant_hits",
            "n_unique_kos",
            "ambiguous_ko_assignment",
            "ribosomal_reference",
        ],
    )

    reference_gene_ids = sorted(hits_by_gene)

    ids_output.parent.mkdir(parents=True, exist_ok=True)
    ids_output.write_text(
        "".join(f"{gene_id}\n" for gene_id in reference_gene_ids),
        encoding="utf-8",
    )

    unique_observed_kos = {hit["ko"] for hit in significant_hits}

    multi_ko_proteins = sum(
        len({hit["ko"] for hit in hits}) > 1 for hits in hits_by_protein.values()
    )

    qc_rows = [
        {
            "sample": sample,
            "n_proteins_total": len(protein_ids),
            "n_reference_kos": len(allowed_kos),
            "n_observed_kos": len(unique_observed_kos),
            "ko_coverage": (
                len(unique_observed_kos) / len(allowed_kos) if allowed_kos else 0.0
            ),
            "n_significant_hits": len(significant_hits),
            "n_ribosomal_proteins": len(hits_by_protein),
            "n_ribosomal_genes": len(hits_by_gene),
            "n_multi_ko_proteins": multi_ko_proteins,
            "n_unmapped_significant_proteins": 0,
        }
    ]

    write_tsv(
        qc_output,
        qc_rows,
        [
            "sample",
            "n_proteins_total",
            "n_reference_kos",
            "n_observed_kos",
            "ko_coverage",
            "n_significant_hits",
            "n_ribosomal_proteins",
            "n_ribosomal_genes",
            "n_multi_ko_proteins",
            "n_unmapped_significant_proteins",
        ],
    )


if __name__ == "__main__":
    main()
