#!/usr/bin/env bash

# Script extracting given in .bed file regions from genomic GCA*.fna files. Downsampling
# Usage:
# ./extract_regions.sh regions.bed genomes extracted_regions

set -euo pipefail

if [[ "$#" -ne 3 ]]; then
    echo "Usage:"
    echo "  $0 <regions.bed> <genome_dir> <out_dir>"
    echo
    echo "Example:"
    echo "  $0 regions.bed genomes extracted_regions"
    exit 1
fi

BED="$1"
GENOME_DIR="$2"
OUT_DIR="$3"

if [[ ! -f "$BED" ]]; then
    echo "ERROR: BED file not found: $BED" >&2
    exit 1
fi

if [[ ! -d "$GENOME_DIR" ]]; then
    echo "ERROR: genome_dir not found: $GENOME_DIR" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

cut -f4 "$BED" | cut -d'|' -f1 | sort -u | while read -r GCA_ID
do
    echo "Processing ${GCA_ID}"

    matches=( "${GENOME_DIR}/${GCA_ID}"*.fna )

    if [[ ! -e "${matches[0]}" ]]; then
        echo "WARNING: genome file not found for ${GCA_ID} in ${GENOME_DIR}" >&2
        continue
    fi

    if [[ "${#matches[@]}" -gt 1 ]]; then
        echo "WARNING: multiple genome files found for ${GCA_ID}, using first:" >&2
        printf '  %s\n' "${matches[@]}" >&2
    fi

    GENOME="${matches[0]}"
    SUB_BED="${OUT_DIR}/${GCA_ID}.bed"
    OUT_FASTA="${OUT_DIR}/${GCA_ID}_regions.fasta"

    if [[ ! -f "$GENOME" ]]; then
        echo "WARNING: genome file not found: $GENOME" >&2
        continue
    fi

    awk -v id="$GCA_ID" -F'\t' '$4 ~ "^"id"\\|" {print}' "$BED" > "$SUB_BED"

    bedtools getfasta \
        -fi "$GENOME" \
        -bed "$SUB_BED" \
        -fo "$OUT_FASTA" \
        -name
done
