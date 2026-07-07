#!/usr/bin/env python3

"""
Extract coding DNA sequences corresponding to ribosomal reference genes.

KofamScan identifies ribosomal proteins using amino-acid sequences, whereas
codon adaptation index and related codon-usage calculations require
nucleotide coding sequences. This script uses the gene identifiers selected
by the KofamScan parsing step to extract the corresponding CDS records from
a genome-specific CDS FASTA file.

Every requested gene identifier must match exactly one FASTA record ID.
The script fails when a selected gene is absent from the CDS FASTA, thereby
preventing silent construction of an incomplete or incorrectly mapped
reference set.

Snakemake inputs
----------------
cds
    FASTA file containing coding DNA sequences for the analysed genome.
gene_ids
    Plain-text file containing one selected ribosomal gene identifier
    per line.

Snakemake outputs
-----------------
cds
    FASTA file containing only coding DNA sequences of the selected
    ribosomal reference genes.

Identifier assumption
---------------------
The first whitespace-delimited token in each CDS FASTA header, exposed as
``SeqRecord.id`` by Biopython, must be identical to the gene identifiers in
``ribosomal_reference_gene_ids.txt``.

If the CDS FASTA uses transcript IDs or CDS IDs instead of gene IDs, an
additional gene-to-CDS mapping must be applied before extraction.

Why this script is required
---------------------------
The functional classification step is performed on protein sequences, but
codon-usage statistics must be calculated from nucleotide CDS sequences.
This script performs the final, validated transition from ribosomal gene
annotations to the nucleotide reference FASTA required by CAI estimation.
"""

from pathlib import Path

from Bio import SeqIO


def main() -> None:
    cds_path = Path(snakemake.input.cds)
    gene_ids_path = Path(snakemake.input.gene_ids)
    output_path = Path(snakemake.output.cds)

    wanted = {
        line.strip()
        for line in gene_ids_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    }

    selected = []
    observed = set()

    for record in SeqIO.parse(cds_path, "fasta"):
        # Zakłada, że pierwszy token nagłówka CDS jest gene_id.
        if record.id in wanted:
            selected.append(record)
            observed.add(record.id)

    missing = wanted - observed

    if missing:
        examples = ", ".join(sorted(missing)[:20])

        raise RuntimeError(
            f"{len(missing)} ribosomal genes were not found in "
            f"{cds_path}. Examples: {examples}"
        )

    if not selected:
        raise RuntimeError("No ribosomal reference CDS were selected")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(selected, output_path, "fasta")


if __name__ == "__main__":
    main()
