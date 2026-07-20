#!/usr/bin/env python3
"""Convert Pfam-A.hmm or Pfam-A.hmm.dat metadata to a compact TSV.

Output columns:
    pfam_id    pfam_description

Accepted records include HMMER fields (ACC/DESC/NAME) and Stockholm general
features (#=GF AC/#=GF DE/#=GF ID). Version suffixes such as PF00001.27 are
removed so the identifiers match annotation tables that use PF00001.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterable, Iterator

PFAM_ID_RE = re.compile(r"\b(PF\d{5})(?:\.\d+)?\b", re.IGNORECASE)


def iter_records(path: Path) -> Iterator[list[str]]:
    record: list[str] = []
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line.strip() == "//":
                if record:
                    yield record
                record = []
            else:
                record.append(line)
    if record:
        yield record


def get_field(record: Iterable[str], *tags: str) -> str | None:
    tag_group = "|".join(re.escape(tag) for tag in tags)
    pattern = re.compile(rf"^(?:#=GF\s+)?(?:{tag_group})\s+(.+?)\s*$")
    for line in record:
        match = pattern.match(line)
        if match:
            return " ".join(match.group(1).split())
    return None


def parse_record(record: list[str]) -> tuple[str, str] | None:
    accession = get_field(record, "ACC", "AC")
    if accession is None:
        return None

    match = PFAM_ID_RE.search(accession)
    if match is None:
        return None

    pfam_id = match.group(1).upper()
    description = get_field(record, "DESC", "DE")
    if not description:
        description = get_field(record, "NAME", "ID")
    if not description:
        return None

    return pfam_id, description


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a PFAM accession-to-description TSV."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    if not args.input.is_file():
        raise FileNotFoundError(f"Input file does not exist: {args.input}")

    descriptions: dict[str, str] = {}
    skipped = 0

    for record in iter_records(args.input):
        parsed = parse_record(record)
        if parsed is None:
            skipped += 1
            continue
        pfam_id, description = parsed
        descriptions.setdefault(pfam_id, description)

    if not descriptions:
        raise RuntimeError(
            "No PFAM records were parsed. Expected ACC/AC and DESC/DE fields."
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("pfam_id", "pfam_description"))
        for pfam_id in sorted(descriptions):
            writer.writerow((pfam_id, descriptions[pfam_id]))

    print(f"Wrote {len(descriptions):,} PFAM descriptions to {args.output}")
    if skipped:
        print(f"Skipped {skipped:,} records without a usable accession/description")


if __name__ == "__main__":
    main()