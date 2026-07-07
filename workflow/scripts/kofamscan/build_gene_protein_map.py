#!/usr/bin/env python3

"""
Build a mapping between protein, gene and CDS identifiers.

The script parses an NCBI-style CDS FASTA file and a corresponding protein
FASTA file. For every CDS record containing a ``[protein_id=...]`` attribute,
it extracts:

- ``protein_id`` — identifier used in the protein FASTA and KofamScan output;
- ``cds_id`` — first token of the CDS FASTA header, equivalent to
  ``Bio.SeqRecord.id``;
- ``gene_id`` — stable gene-level identifier selected using the following
  priority:

    1. ``locus_tag``,
    2. ``gene``,
    3. ``protein_id`` as a fallback.

The resulting table allows downstream workflow steps to connect KofamScan
protein-level annotations with gene-level analysis tables and the original
nucleotide CDS records.

The script validates that:

1. protein identifiers in the protein FASTA are unique;
2. protein identifiers in the CDS FASTA are unique;
3. each mapped protein exists in the corresponding protein FASTA when
   strict matching is enabled;
4. one protein identifier does not map to multiple CDS or genes.

Snakemake inputs
----------------
cds
    NCBI-style nucleotide CDS FASTA file containing ``protein_id`` and
    preferably ``locus_tag`` annotations in record descriptions.
proteins
    Protein FASTA file used as input for KofamScan.

Snakemake output
----------------
mapping
    Tab-separated mapping table containing sample, gene ID, CDS ID,
    protein ID and identifier provenance.

Snakemake parameters
--------------------
sample
    Sample identifier written to every output row.
strict_proteome_match
    When true, fail if a CDS protein identifier is absent from the protein
    FASTA. When false, retain the mapping and mark the protein as absent.

Why this script is required
---------------------------
KofamScan reports annotations using protein identifiers, while codon-usage
analysis uses nucleotide CDS records and canonical tables operate at the
gene level. These identifiers are different in NCBI FASTA files and cannot
be safely treated as interchangeable.
"""

import csv
import gzip
import re
from pathlib import Path

from Bio import SeqIO

ATTRIBUTE_PATTERN = re.compile(r"\[([^=\]]+)=([^\]]*)\]")


def open_fasta(path: Path):
    """Open an uncompressed or gzip-compressed FASTA file in text mode."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_attributes(description: str) -> dict[str, str]:
    """Extract NCBI-style ``[key=value]`` attributes from a FASTA header."""
    return {
        key.strip(): value.strip()
        for key, value in ATTRIBUTE_PATTERN.findall(description)
    }


def clean_text(value: str) -> str:
    """Remove tab and newline characters before writing a TSV field."""
    return value.replace("\t", " ").replace("\n", " ").strip()


def read_protein_ids(path: Path) -> set[str]:
    """Read and validate identifiers from the protein FASTA."""
    protein_ids: set[str] = set()
    duplicates: set[str] = set()

    with open_fasta(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in protein_ids:
                duplicates.add(record.id)
            protein_ids.add(record.id)

    if duplicates:
        preview = ", ".join(sorted(duplicates)[:20])
        raise ValueError("Duplicate identifiers in protein FASTA: " f"{preview}")

    if not protein_ids:
        raise ValueError(f"No protein records found in {path}")

    return protein_ids


def select_gene_id(
    attributes: dict[str, str],
    protein_id: str,
) -> tuple[str, str]:
    """Select a gene identifier and report its source field."""
    locus_tag = attributes.get("locus_tag", "").strip()
    if locus_tag:
        return locus_tag, "locus_tag"

    gene_name = attributes.get("gene", "").strip()
    if gene_name:
        return gene_name, "gene"

    return protein_id, "protein_id_fallback"


def main() -> None:
    cds_path = Path(snakemake.input.cds)
    proteins_path = Path(snakemake.input.proteins)
    output_path = Path(snakemake.output.mapping)

    sample = str(snakemake.params.sample)
    strict_match = bool(
        getattr(
            snakemake.params,
            "strict_proteome_match",
            True,
        )
    )

    protein_fasta_ids = read_protein_ids(proteins_path)

    rows: list[dict[str, object]] = []
    seen_protein_ids: dict[str, tuple[str, str]] = {}
    missing_from_proteome: list[str] = []
    records_without_protein_id: list[str] = []

    with open_fasta(cds_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            attributes = parse_attributes(record.description)

            protein_id = attributes.get("protein_id", "").strip()

            if not protein_id:
                records_without_protein_id.append(record.id)
                continue

            gene_id, gene_id_source = select_gene_id(
                attributes,
                protein_id,
            )

            previous = seen_protein_ids.get(protein_id)
            current = (gene_id, record.id)

            if previous is not None and previous != current:
                raise ValueError(
                    f"Protein {protein_id} maps to multiple records: "
                    f"{previous} and {current}"
                )

            seen_protein_ids[protein_id] = current

            present_in_proteome = protein_id in protein_fasta_ids

            if not present_in_proteome:
                missing_from_proteome.append(protein_id)

            rows.append(
                {
                    "sample": sample,
                    "gene_id": gene_id,
                    "gene_id_source": gene_id_source,
                    "locus_tag": attributes.get("locus_tag", ""),
                    "gene_name": attributes.get("gene", ""),
                    "protein_id": protein_id,
                    "cds_id": record.id,
                    "present_in_proteome": str(present_in_proteome).lower(),
                    "cds_description": clean_text(record.description),
                }
            )

    if not rows:
        raise RuntimeError(
            f"No CDS-to-protein mappings could be created from {cds_path}"
        )

    if strict_match and missing_from_proteome:
        preview = ", ".join(sorted(set(missing_from_proteome))[:20])

        raise RuntimeError(
            f"{len(set(missing_from_proteome))} CDS protein identifiers "
            f"are absent from {proteins_path}. Examples: {preview}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")

    fieldnames = [
        "sample",
        "gene_id",
        "gene_id_source",
        "locus_tag",
        "gene_name",
        "protein_id",
        "cds_id",
        "present_in_proteome",
        "cds_description",
    ]

    with temporary_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    temporary_path.replace(output_path)

    print(f"sample: {sample}")
    print(f"protein_fasta_records: {len(protein_fasta_ids)}")
    print(f"mapped_cds_records: {len(rows)}")
    print("cds_without_protein_id: " f"{len(records_without_protein_id)}")
    print("cds_proteins_missing_from_proteome: " f"{len(set(missing_from_proteome))}")


if __name__ == "__main__":
    main()
