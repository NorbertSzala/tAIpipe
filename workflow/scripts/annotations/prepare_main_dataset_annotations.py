#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv
import re
import sys

GO_BRACE_RE = re.compile(r"\{(\d{7})\}")
PFAM_RE = re.compile(r"PF\d{5}(?:\.\d+)?")

COLUMNS = [
    "accession",
    "protein_id",
    "protein_length_main_dataset",
    "pfam_present",
    "pfam_domain_bins_10",
    "tm_count",
    "tm_total_length",
    "tm_bins_10",
    "signal_peptide_present",
    "lcr_count",
    "lcr_total_length",
    "lcr_bins_3",
    "lcr_bins_10",
    "lcr_sequences_n",
    "lcr_sequences_middle",
    "lcr_sequences_c",
    "pfam_lcr_overlap_terms",
    "pfam_terms",
    "go_terms_raw",
]


def as_bool01(x):
    x = str(x).strip()
    if x in {"1", "true", "TRUE", "True", "yes", "Y"}:
        return "TRUE"
    if x in {"0", "false", "FALSE", "False", "no", "N", "", "NA"}:
        return "FALSE"
    return "NA"


def normalize_pfam(x):
    hits = PFAM_RE.findall(str(x))
    return ";".join(hits) if hits else "NA"


def normalize_go(x):
    ids = sorted(set(f"GO:{m}" for m in GO_BRACE_RE.findall(str(x))))
    return ";".join(ids) if ids else "NA"


def load_sample_accessions(samples_tsv):
    mapping = {}
    with open(samples_tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"sample", "accession", "include"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"samples.tsv missing columns: {sorted(missing)}")

        for row in reader:
            include = str(row["include"]).strip().lower() in {
                "true",
                "1",
                "yes",
                "y",
                "t",
            }
            if not include:
                continue
            acc = row["accession"].strip()
            sample = row["sample"].strip()
            if acc:
                mapping[acc] = sample

    return mapping


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--main-dataset", required=True)
    p.add_argument("--samples", default="config/samples.tsv")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    accession_to_sample = load_sample_accessions(args.samples)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    n_in = 0
    n_out = 0
    n_missing_sample = 0

    with open(args.main_dataset, newline="") as fin, out.open("w", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.DictWriter(
            fout,
            delimiter="\t",
            fieldnames=[
                "sample",
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
                "pfam_lcr_overlap_terms",
                "go_terms",
                "annotation_source",
            ],
        )
        writer.writeheader()

        for row in reader:
            n_in += 1

            if len(row) < len(COLUMNS):
                row = row + [""] * (len(COLUMNS) - len(row))

            rec = dict(zip(COLUMNS, row[: len(COLUMNS)]))
            accession = rec["accession"].strip()
            sample = accession_to_sample.get(accession)

            if not sample:
                n_missing_sample += 1
                continue

            protein_id = rec["protein_id"].strip()
            if not protein_id:
                continue

            tm_count = int(rec["tm_count"]) if rec["tm_count"].strip().isdigit() else 0
            lcr_count = (
                int(rec["lcr_count"]) if rec["lcr_count"].strip().isdigit() else 0
            )

            pfam_terms = normalize_pfam(rec["pfam_terms"])
            pfam_lcr_overlap_terms = normalize_pfam(rec["pfam_lcr_overlap_terms"])
            go_terms = normalize_go(rec["go_terms_raw"])

            writer.writerow(
                {
                    "sample": sample,
                    "protein_id": protein_id,
                    "protein_length_main_dataset": rec["protein_length_main_dataset"],
                    "signal_peptide_present": as_bool01(rec["signal_peptide_present"]),
                    "tm_present": "TRUE" if tm_count > 0 else "FALSE",
                    "tm_count": tm_count,
                    "tm_total_length": rec["tm_total_length"] or "0",
                    "lcr_present": "TRUE" if lcr_count > 0 else "FALSE",
                    "lcr_count": lcr_count,
                    "lcr_total_length": rec["lcr_total_length"] or "0",
                    "pfam_present": as_bool01(rec["pfam_present"]),
                    "pfam_terms": pfam_terms,
                    "pfam_lcr_overlap_terms": pfam_lcr_overlap_terms,
                    "go_terms": go_terms,
                    "annotation_source": "Main_Dataset.tsv",
                }
            )
            n_out += 1

    print(f"input_rows={n_in}", file=sys.stderr)
    print(f"output_rows={n_out}", file=sys.stderr)
    print(f"rows_skipped_missing_sample={n_missing_sample}", file=sys.stderr)


if __name__ == "__main__":
    main()
