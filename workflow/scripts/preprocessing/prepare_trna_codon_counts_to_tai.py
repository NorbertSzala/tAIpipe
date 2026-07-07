#!/usr/bin/env python

"""
# Converts cleaned tRNA predictions into counts grouped by amino acid and anticodon. It filters initiator, undetermined, and malformed records and writes the normalized tRNA copy-number profile expected by the tAI calculation script.

Input:
    cleaned tRNAscan-SE output .tsv

Expected columns:
    seq_name
    trna_number
    begin
    end
    trna_type
    anticodon
    intron_begin
    intron_end
    score
    note
    strand
    start
    stop
    has_intron
    is_pseudo

Output:
    TSV with tRNA gene copy numbers in format:
        anticodon_id    count
        Ala-CGC         5

Usage:
    ./scripts/prepare_trna_counts_to_tai.py \
        -I results/trnascan_cleaned.tsv \
        -O results/trna_levels.tsv
"""

# ----------------------
# --- Import modules ---
# ----------------------

from pathlib import Path
import argparse
import pandas as pd


# -----------------
# --- Arguments ---
# -----------------
def arguments():
    p = argparse.ArgumentParser(
        prog="prepare_trna_codon_counts_to_tai",
        description="Script to extract Aminoacid-Anticodon counts",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        conflict_handler="error",
        add_help=True,
    )

    p.add_argument(
        "-I",
        "--input",
        help="Path to input file (cleaned tRNAscan-SE processed by convert_trnascance_output_to_tsv.py)",
        required=True,
        type=Path,
    )

    p.add_argument(
        "-O",
        "--output",
        help="Path to ouptut .tsv file",
        required=True,
        type=Path,
    )

    return p.parse_args()


# ------------------------
# --- Helper functions ---
# ------------------------
def read_cleaned_trnascan(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=0)
    required = {"trna_type", "anticodon"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}",
            f"Available columns: {list(df.columns)}",
        )

    return df


def count_aminoacid(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["aminoacid-anticodon"] = df["trna_type"] + "-" + df["anticodon"]

    counts = (
        df.groupby("aminoacid-anticodon", as_index=False)
        .size()
        .rename(columns={"size": "count"})
        .sort_values(["aminoacid-anticodon"])
    )
    return counts


def write_counts(counts: pd.DataFrame, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output, sep="\t", index=False)


# ---------------------
# --- Main function ---
# ---------------------
def main():
    args = arguments()
    df = read_cleaned_trnascan(args.input)

    if "is_pseudo" in df.columns:
        df = df[df["is_pseudo"].astype(str).str.lower() == "false"].copy()

    trna_type = df["trna_type"].astype(str).str.strip().str.lower()
    anticodon = (
        df["anticodon"]
        .astype(str)
        .str.strip()
        .str.upper()
        .str.replace("U", "T", regex=False)
    )

    df = df[
        ~trna_type.isin({"imet", "fmet", "und"})
        & anticodon.str.fullmatch(r"[ACGT]{3}", na=False)
    ].copy()

    df["anticodon"] = anticodon.loc[df.index]

    counts = count_aminoacid(df)
    write_counts(counts, args.output)

    print(f"input_rows_after_filtering: {len(df)}")
    print(f"unique_aminoacid_anticodon: {len(counts)}")
    print(f"output: {args.output}")


if __name__ == "__main__":
    main()
