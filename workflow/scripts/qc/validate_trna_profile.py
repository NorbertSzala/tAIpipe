#!/usr/bin/env python3
"""
# Evaluates the structural validity and biological completeness of a sample's tRNA copy-number profile. It reports total and elongator tRNA counts, anticodon diversity, amino-acid coverage, invalid identifiers, duplicate records, and an overall PASS, WARN, or FAIL status based on configurable thresholds.

Validate a tRNA gene-copy-number profile before tAI calculation.

Expected input columns:
    aminoacid-anticodon    count

The first column may also have another name; the first two columns are treated
as the tRNA identifier and its count. The script always writes a one-row QC TSV.
Execution behaviour is controlled by --mode:
    error  - return exit code 1 for WARN or FAIL
    warn   - return exit code 0 and print non-PASS status to stderr
    ignore - return exit code 0 without an additional warning
"""

# ----------------------
# --- Import modules ---
# ----------------------
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import pandas as pd

# -----------------------
# --- Constant values ---
# -----------------------
STANDARD_AMINO_ACIDS = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
INITIATOR_TYPES = {"IMET", "FMET", "INITIATOR", "INITIATORMET"}
UNDETERMINED_TYPES = {"UND", "UNKNOWN", "NNN"}
IDENTIFIER_PATTERN = re.compile(r"^(?P<aa>[^-]+)-(?P<anticodon>[ACGTU]{3}|NNN)$", re.I)


# -----------------
# --- Arguments ---
# -----------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate a tRNA copy-number profile used by cubar/tAI.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-I", "--input", required=True, type=Path)
    parser.add_argument("-O", "--output", required=True, type=Path)
    parser.add_argument("--sample", default="")
    parser.add_argument("--mode", choices=("error", "warn", "ignore"), default="warn")
    parser.add_argument("--min-total-trnas", type=int, default=20)
    parser.add_argument("--min-unique-anticodons", type=int, default=15)
    parser.add_argument("--min-amino-acids", type=int, default=18)
    return parser.parse_args()


# ------------------------
# --- Helper functions ---
# ------------------------
def infer_sample(path: Path) -> str:
    suffixes = (
        "_aaa_counts.tsv",
        "_trna_profile.tsv",
        ".tsv",
    )
    name = path.name
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def join_reasons(reasons: list[str]) -> str:
    return ";".join(dict.fromkeys(reasons)) if reasons else ""


def write_qc(output: Path, row: dict[str, object]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(output, sep="\t", index=False, na_rep="NA")


# ---------------------
# --- Main function ---
# ---------------------
def main() -> int:
    args = parse_args()
    sample = args.sample.strip() or infer_sample(args.input)

    if not args.input.exists():
        raise FileNotFoundError(f"Input tRNA count file does not exist: {args.input}")
    if args.input.stat().st_size == 0:
        raise ValueError(f"Input tRNA count file is empty: {args.input}")

    df = pd.read_csv(args.input, sep="\t", dtype=str)
    if df.shape[1] < 2:
        raise ValueError(
            "tRNA count file must contain at least two columns: "
            "aminoacid-anticodon and count"
        )

    id_column, count_column = df.columns[:2]
    data = df[[id_column, count_column]].copy()
    data.columns = ["trna_id", "count_raw"]
    data["trna_id"] = data["trna_id"].fillna("").str.strip()
    data["count"] = pd.to_numeric(data["count_raw"], errors="coerce")

    fail_reasons: list[str] = []
    warn_reasons: list[str] = []

    n_rows = len(data)
    n_missing_ids = int(data["trna_id"].eq("").sum())
    n_non_numeric_counts = int(data["count"].isna().sum())
    n_negative_counts = int((data["count"].fillna(0) < 0).sum())
    n_non_integer_counts = int(
        data["count"]
        .dropna()
        .map(lambda value: not math.isclose(value, round(value)))
        .sum()
    )
    n_duplicate_ids = int(data["trna_id"].duplicated(keep=False).sum())

    if n_rows == 0:
        fail_reasons.append("empty_profile")
    if n_missing_ids:
        fail_reasons.append("missing_trna_identifiers")
    if n_non_numeric_counts:
        fail_reasons.append("non_numeric_counts")
    if n_negative_counts:
        fail_reasons.append("negative_counts")
    if n_non_integer_counts:
        fail_reasons.append("non_integer_counts")
    if n_duplicate_ids:
        fail_reasons.append("duplicated_trna_identifiers")

    parsed = data["trna_id"].str.extract(IDENTIFIER_PATTERN)
    data["amino_acid"] = parsed["aa"].fillna("").str.upper().str.strip()
    data["anticodon"] = (
        parsed["anticodon"].fillna("").str.upper().str.replace("U", "T", regex=False)
    )
    n_invalid_identifiers = int(parsed["aa"].isna().sum())
    if n_invalid_identifiers:
        fail_reasons.append("invalid_aminoacid_anticodon_identifiers")

    valid_numeric = data[
        data["count"].notna()
        & (data["count"] >= 0)
        & data["amino_acid"].ne("")
        & data["anticodon"].ne("")
    ].copy()

    valid_numeric["count"] = valid_numeric["count"].round().astype(int)
    valid_numeric["is_initiator"] = valid_numeric["amino_acid"].isin(INITIATOR_TYPES)
    valid_numeric["is_undetermined"] = valid_numeric["amino_acid"].isin(
        UNDETERMINED_TYPES
    ) | valid_numeric["anticodon"].eq("NNN")
    valid_numeric["is_standard_amino_acid"] = valid_numeric["amino_acid"].isin(
        STANDARD_AMINO_ACIDS
    )

    elongator = valid_numeric[
        ~valid_numeric["is_initiator"]
        & ~valid_numeric["is_undetermined"]
        & valid_numeric["is_standard_amino_acid"]
        & (valid_numeric["count"] > 0)
    ].copy()

    n_total_trnas = int(valid_numeric["count"].sum()) if not valid_numeric.empty else 0
    n_elongator_trnas = int(elongator["count"].sum()) if not elongator.empty else 0
    n_initiator_trnas = (
        int(valid_numeric.loc[valid_numeric["is_initiator"], "count"].sum())
        if not valid_numeric.empty
        else 0
    )
    n_undetermined_trnas = (
        int(valid_numeric.loc[valid_numeric["is_undetermined"], "count"].sum())
        if not valid_numeric.empty
        else 0
    )
    n_unique_aa_anticodon = int(elongator["trna_id"].nunique())
    n_unique_anticodons = int(elongator["anticodon"].nunique())
    n_amino_acids = int(elongator["amino_acid"].nunique())
    missing_amino_acids = sorted(
        STANDARD_AMINO_ACIDS - set(elongator["amino_acid"].unique())
    )

    if n_elongator_trnas == 0:
        fail_reasons.append("no_elongator_trnas")

    if n_total_trnas < args.min_total_trnas:
        warn_reasons.append(f"total_trnas_below_{args.min_total_trnas}")
    if n_unique_anticodons < args.min_unique_anticodons:
        warn_reasons.append(f"unique_anticodons_below_{args.min_unique_anticodons}")
    if n_amino_acids < args.min_amino_acids:
        warn_reasons.append(f"amino_acids_below_{args.min_amino_acids}")
    if n_initiator_trnas > 0:
        warn_reasons.append("initiator_trnas_present_in_input")
    if n_undetermined_trnas > 0:
        warn_reasons.append("undetermined_trnas_present_in_input")

    if fail_reasons:
        qc_status = "FAIL"
    elif warn_reasons:
        qc_status = "WARN"
    else:
        qc_status = "PASS"

    all_reasons = fail_reasons + warn_reasons
    row = {
        "sample": sample,
        "input_file": str(args.input),
        "n_rows": n_rows,
        "n_total_trnas": n_total_trnas,
        "n_elongator_trnas": n_elongator_trnas,
        "n_initiator_trnas": n_initiator_trnas,
        "n_undetermined_trnas": n_undetermined_trnas,
        "n_unique_aa_anticodon": n_unique_aa_anticodon,
        "n_unique_anticodons": n_unique_anticodons,
        "n_amino_acids": n_amino_acids,
        "missing_amino_acids": ";".join(missing_amino_acids),
        "n_missing_ids": n_missing_ids,
        "n_invalid_identifiers": n_invalid_identifiers,
        "n_non_numeric_counts": n_non_numeric_counts,
        "n_negative_counts": n_negative_counts,
        "n_non_integer_counts": n_non_integer_counts,
        "n_duplicate_ids": n_duplicate_ids,
        "min_total_trnas": args.min_total_trnas,
        "min_unique_anticodons": args.min_unique_anticodons,
        "min_amino_acids": args.min_amino_acids,
        "qc_status": qc_status,
        "qc_reasons": join_reasons(all_reasons),
    }
    write_qc(args.output, row)

    print(f"sample: {sample}")
    print(f"qc_status: {qc_status}")
    print(f"n_total_trnas: {n_total_trnas}")
    print(f"n_elongator_trnas: {n_elongator_trnas}")
    print(f"n_unique_anticodons: {n_unique_anticodons}")
    print(f"n_amino_acids: {n_amino_acids}")
    print(f"qc_reasons: {join_reasons(all_reasons) or 'none'}")
    print(f"output: {args.output}")

    if qc_status != "PASS":
        message = (
            f"tRNA profile QC for {sample}: {qc_status}: {join_reasons(all_reasons)}"
        )
        if args.mode == "warn":
            print(f"WARNING: {message}", file=sys.stderr)
        elif args.mode == "error":
            print(f"ERROR: {message}", file=sys.stderr)
            return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
