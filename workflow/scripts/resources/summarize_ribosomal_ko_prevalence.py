#!/usr/bin/env python3
"""Summarise prevalence and copy-number consistency of ribosomal KOfam markers."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", default="config/samples.tsv")
    ap.add_argument("--results-root", default="results/per_genome")
    ap.add_argument(
        "--output",
        default="results/tables/ribosomal_ko_prevalence.tsv",
    )
    args = ap.parse_args()

    with Path(args.samples).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    expected = {
        row["sample"]
        for row in rows
        if "include" not in row or truthy(row["include"])
    }
    if not expected:
        raise SystemExit("No included samples found in the sample sheet.")

    # hits[KO][sample] = set of distinct protein IDs assigned to that KO.
    hits: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    definitions: dict[str, Counter[str]] = defaultdict(Counter)
    missing_files: list[str] = []

    for sample in sorted(expected):
        path = (
            Path(args.results_root)
            / sample
            / "kofamscan"
            / "ribosome_significant_hits.tsv"
        )
        if not path.exists():
            missing_files.append(sample)
            continue

        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                ko = row.get("ko", "").strip()
                protein = row.get("protein_id", "").strip()
                definition = row.get("definition", "").strip()
                if ko and protein:
                    hits[ko][sample].add(protein)
                    if definition:
                        definitions[ko][definition] += 1

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "ko",
        "definition",
        "n_expected_genomes",
        "n_genomes_present",
        "prevalence",
        "n_genomes_exactly_one_hit",
        "single_hit_prevalence",
        "n_genomes_multiple_hits",
        "passes_80pct_presence",
        "passes_90pct_presence",
        "passes_80pct_single_hit",
        "passes_90pct_single_hit",
    ]

    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()

        for ko in sorted(hits):
            per_sample = hits[ko]
            present = len(per_sample)
            single = sum(len(ids) == 1 for ids in per_sample.values())
            multiple = sum(len(ids) > 1 for ids in per_sample.values())
            n = len(expected)
            definition = (
                definitions[ko].most_common(1)[0][0]
                if definitions[ko]
                else ""
            )
            writer.writerow(
                {
                    "ko": ko,
                    "definition": definition,
                    "n_expected_genomes": n,
                    "n_genomes_present": present,
                    "prevalence": f"{present / n:.6f}",
                    "n_genomes_exactly_one_hit": single,
                    "single_hit_prevalence": f"{single / n:.6f}",
                    "n_genomes_multiple_hits": multiple,
                    "passes_80pct_presence": present / n >= 0.80,
                    "passes_90pct_presence": present / n >= 0.90,
                    "passes_80pct_single_hit": single / n >= 0.80,
                    "passes_90pct_single_hit": single / n >= 0.90,
                }
            )

    print(f"Expected genomes: {len(expected)}")
    print(f"Genomes with missing hit table: {len(missing_files)}")
    print(f"Distinct KOs observed: {len(hits)}")
    print(
        "KOs present in >=80% / >=90%: "
        f"{sum(len(v) / len(expected) >= 0.80 for v in hits.values())} / "
        f"{sum(len(v) / len(expected) >= 0.90 for v in hits.values())}"
    )
    print(
        "KOs with exactly one hit in >=80% / >=90%: "
        f"{sum(sum(len(x) == 1 for x in v.values()) / len(expected) >= 0.80 for v in hits.values())} / "
        f"{sum(sum(len(x) == 1 for x in v.values()) / len(expected) >= 0.90 for v in hits.values())}"
    )
    print(f"Wrote: {out}")

    if missing_files:
        missing_out = out.with_name(out.stem + ".missing_samples.txt")
        missing_out.write_text("\n".join(missing_files) + "\n", encoding="utf-8")
        print(f"Missing-sample list: {missing_out}")


if __name__ == "__main__":
    main()
