#!/usr/bin/env python3
"""
Convert the raw Main_Dataset.tsv annotation table into a pipeline-ready TSV.

Input assumptions
-----------------
The raw file has no header and contains the columns described in the project
notes: assembly accession, protein accession, protein length, PFAM/LCR/TM/signal
features, PFAM domains overlapping LCRs, all PFAM domains, and GO annotations
created through PFAM/PFAM2GO-like mapping.

Why this script exists
----------------------
build_gene_features.R expects an annotation table keyed by sample + protein_id
or sample + seq_id. The raw Main_Dataset.tsv is keyed by NCBI assembly accession
and protein accession. This script resolves assembly accession -> sample using
config/samples.tsv and normalizes feature columns into stable machine-readable
fields used by plotting/statistical scripts.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Iterable

GO_BRACE_RE = re.compile(r"\{(\d{7})\}")
GO_ID_RE = re.compile(r"GO:\d{7}", re.IGNORECASE)
PFAM_RE = re.compile(r"PF\d{5}(?:\.\d+)?")

RAW_COLUMNS = [
    "accession",
    "protein_id",
    "protein_length_main_dataset",
    "pfam_present_raw",
    "pfam_domain_bins_10",
    "tm_count",
    "tm_total_length",
    "tm_bins_10",
    "signal_peptide_present_raw",
    "lcr_count",
    "lcr_total_length",
    "lcr_bins_3",
    "lcr_bins_10",
    "lcr_sequences_n",
    "lcr_sequences_middle",
    "lcr_sequences_c",
    "pfam_lcr_overlap_terms_raw",
    "pfam_terms_raw",
    "go_terms_raw",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare extended Main_Dataset annotations for tAIpipe."
    )
    parser.add_argument("--main-dataset", required=True, help="Raw Main_Dataset.tsv")
    parser.add_argument("--samples", default="config/samples.tsv", help="samples.tsv")
    parser.add_argument("--output", required=True, help="Output annotation TSV")
    parser.add_argument(
        "--unmapped-output",
        default=None,
        help="Optional TSV with rows skipped because assembly accession was not in samples.tsv",
    )
    return parser.parse_args()


def clean_text(value: object) -> str:
    if value is None:
        return ""
    return str(value).strip()


def as_int(value: object, default: int = 0) -> int:
    text = clean_text(value)
    try:
        return int(float(text))
    except Exception:
        return default


def as_float_text(value: object, default: str = "0") -> str:
    text = clean_text(value)
    if text == "":
        return default
    try:
        return str(float(text)).rstrip("0").rstrip(".")
    except Exception:
        return default


def bool01(value: object) -> str:
    text = clean_text(value).lower()
    if text in {"1", "true", "t", "yes", "y"}:
        return "TRUE"
    if text in {"0", "false", "f", "no", "n", "", "na", "nan"}:
        return "FALSE"
    return "NA"


def unique_sorted(values: Iterable[str]) -> list[str]:
    cleaned = sorted({v.strip() for v in values if v and v.strip()})
    return cleaned


def normalize_pfam(value: object) -> str:
    hits = unique_sorted(PFAM_RE.findall(clean_text(value)))
    return ";".join(hits) if hits else "NA"


def count_semicolon(value: str) -> int:
    if not value or value == "NA":
        return 0
    return len([x for x in value.split(";") if x])


def normalize_go(value: object) -> str:
    text = clean_text(value)
    brace_hits = [f"GO:{m}" for m in GO_BRACE_RE.findall(text)]
    explicit_hits = [m.upper() for m in GO_ID_RE.findall(text)]
    hits = unique_sorted(brace_hits + explicit_hits)
    return ";".join(hits) if hits else "NA"


def load_accession_to_sample(samples_path: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with open(samples_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "accession", "include"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"samples.tsv missing required columns: {sorted(missing)}")

        for row in reader:
            include = clean_text(row.get("include", "")).lower() in {
                "true", "1", "yes", "y", "t"
            }
            if not include:
                continue
            accession = clean_text(row.get("accession"))
            sample = clean_text(row.get("sample"))
            if accession and sample:
                mapping[accession] = sample
    return mapping


def main() -> None:
    args = parse_args()
    accession_to_sample = load_accession_to_sample(args.samples)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    unmapped_rows = []
    n_in = n_out = n_no_sample = n_no_protein = 0

    out_fields = [
        "sample",
        "accession",
        "protein_id",
        "protein_length_main_dataset",
        "signal_peptide_present",
        "tm_present",
        "tm_count",
        "tm_total_length",
        "lcr_present",
        "lcr_count",
        "lcr_total_length",
        "pfam_present",
        "pfam_terms",
        "pfam_count",
        "pfam_lcr_overlap_terms",
        "pfam_lcr_overlap_count",
        "pfam_lcr_overlap_present",
        "go_terms",
        "go_term_count",
        "annotation_source",
    ]

    with open(args.main_dataset, newline="") as fin, output.open("w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, delimiter="\t", fieldnames=out_fields)
        writer.writeheader()

        for row in reader:
            n_in += 1
            if len(row) < len(RAW_COLUMNS):
                row = row + [""] * (len(RAW_COLUMNS) - len(row))
            rec = dict(zip(RAW_COLUMNS, row[: len(RAW_COLUMNS)]))

            accession = clean_text(rec["accession"])
            protein_id = clean_text(rec["protein_id"])
            sample = accession_to_sample.get(accession)

            if not sample:
                n_no_sample += 1
                unmapped_rows.append(row)
                continue
            if not protein_id:
                n_no_protein += 1
                continue

            tm_count = as_int(rec["tm_count"])
            lcr_count = as_int(rec["lcr_count"])
            pfam_terms = normalize_pfam(rec["pfam_terms_raw"])
            pfam_lcr_terms = normalize_pfam(rec["pfam_lcr_overlap_terms_raw"])
            go_terms = normalize_go(rec["go_terms_raw"])

            writer.writerow(
                {
                    "sample": sample,
                    "accession": accession,
                    "protein_id": protein_id,
                    "protein_length_main_dataset": as_int(
                        rec["protein_length_main_dataset"], default=0
                    ),
                    "signal_peptide_present": bool01(rec["signal_peptide_present_raw"]),
                    "tm_present": "TRUE" if tm_count > 0 else "FALSE",
                    "tm_count": tm_count,
                    "tm_total_length": as_float_text(rec["tm_total_length"]),
                    "lcr_present": "TRUE" if lcr_count > 0 else "FALSE",
                    "lcr_count": lcr_count,
                    "lcr_total_length": as_float_text(rec["lcr_total_length"]),
                    "pfam_present": bool01(rec["pfam_present_raw"]),
                    "pfam_terms": pfam_terms,
                    "pfam_count": count_semicolon(pfam_terms),
                    "pfam_lcr_overlap_terms": pfam_lcr_terms,
                    "pfam_lcr_overlap_count": count_semicolon(pfam_lcr_terms),
                    "pfam_lcr_overlap_present": (
                        "TRUE" if count_semicolon(pfam_lcr_terms) > 0 else "FALSE"
                    ),
                    "go_terms": go_terms,
                    "go_term_count": count_semicolon(go_terms),
                    "annotation_source": "Main_Dataset.tsv",
                }
            )
            n_out += 1

    if args.unmapped_output:
        unmapped = Path(args.unmapped_output)
        unmapped.parent.mkdir(parents=True, exist_ok=True)
        with unmapped.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(RAW_COLUMNS)
            writer.writerows(unmapped_rows)

    print(f"input_rows={n_in}", file=sys.stderr)
    print(f"output_rows={n_out}", file=sys.stderr)
    print(f"rows_skipped_missing_sample={n_no_sample}", file=sys.stderr)
    print(f"rows_skipped_missing_protein={n_no_protein}", file=sys.stderr)


if __name__ == "__main__":
    main()
