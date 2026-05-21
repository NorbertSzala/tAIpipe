#!/usr/bin/env bash

# Script extracting first contig from given genome. Downsampling (choosing 10k sequences) was inaccurate for testing tRNAscan-SE

set -euo pipefail

if [[ "$#" -ne 2 ]]; then
    echo "Usage:"
    echo "  $0 <genome_dir> <out_dir>"
    echo
    echo "Example:"
    echo "  $0 data/raw/genome data/downsampled/first_contigs"
    exit 1
fi

GENOME_DIR="$1"
OUT_DIR="$2"

mkdir -p "$OUT_DIR"


for genome in "$GENOME_DIR"/*.fna
do
    filename=$(basename "$genome")

    accession=$(echo "$filename" | grep -oE '^GCA_[0-9]+\.[0-9]+')

    if [[ -z "$accession" ]]; then
        echo "WARNING: could not extract accession from: $filename" >&2
        continue
    fi

    out="${OUT_DIR}/${accession}_first_contig.fna"

    echo "Processing $genome -> $out"

    awk '
        /^>/ {
            if (seen == 1) {
                exit
            }
            seen = 1
        }
        seen == 1 {
            print
        }
    ' "$genome" > "$out"
done