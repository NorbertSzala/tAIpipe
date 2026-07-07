#!/usr/bin/env python3

"""
Prepare a restricted KofamScan database containing cytosolic
eukaryotic ribosomal protein profiles.

The script reads a locally stored KEGG BRITE ko03011 hierarchy and extracts
KEGG Orthology identifiers from the following branches:

    Ribosomal proteins
      └── Eukaryotes
          ├── Large subunit
          └── Small subunit

For every selected KO identifier, the script verifies that:

1. the KO is present in the local KOfam ``ko_list`` file;
2. the corresponding ``Kxxxxx.hmm`` profile exists;
3. a predefined KOfam score threshold is available, when
   ``require_thresholds`` is enabled.

The script produces two files:

``cytosolic_eukaryotic_ribosome_kos.tsv``
    A human-readable table containing the selected KO identifiers,
    ribosomal subunit assignments, BRITE annotations, KOfam thresholds,
    score types, profile types and functional definitions.

``cytosolic_eukaryotic_ribosome.hal``
    A KofamScan profile-list file containing relative paths to the selected
    HMM profiles. This file can be passed directly to ``exec_annotation``
    through the ``--profile`` option.

Snakemake inputs
----------------
brite_json
    Local KEGG BRITE ko03011 JSON file restricted to the cytosolic
    eukaryotic large and small ribosomal subunits.
ko_list
    KOfam metadata table containing KO-specific thresholds and definitions.

Snakemake outputs
-----------------
ko_table
    Tab-separated table describing the selected ribosomal KO profiles.
hal
    KofamScan HAL file listing the corresponding HMM profiles.

Snakemake parameters
--------------------
profiles_dir
    Directory containing individual KOfam ``Kxxxxx.hmm`` profiles.
require_thresholds
    Whether the script should fail when a selected KO lacks a predefined
    KOfam score threshold.

Why this script is required
---------------------------
KofamScan normally searches thousands of KO profiles. The downstream
analysis requires only cytosolic eukaryotic ribosomal proteins. Restricting
the profile set reduces runtime, prevents mitochondrial and prokaryotic
ribosomal proteins from entering the reference set, and creates an explicit,
inspectable definition of the ribosomal KO reference.
"""

import csv
import json
import os
import re
from pathlib import Path
from typing import Any, Iterator

KO_PATTERN = re.compile(r"\bK\d{5}\b")


def clean_label(value: str) -> str:
    """Remove optional numeric prefix from a BRITE node label."""
    return re.sub(r"^\d+\s+", "", value).strip()


def walk_tree(
    node: dict[str, Any],
    path: tuple[str, ...] = (),
) -> Iterator[tuple[dict[str, Any], tuple[str, ...]]]:
    name = str(node.get("name", "")).strip()
    current_path = path + (name,)

    yield node, current_path

    for child in node.get("children", []):
        yield from walk_tree(child, current_path)


def find_nodes(
    root: dict[str, Any],
    expected_label: str,
) -> list[dict[str, Any]]:
    expected = expected_label.casefold()

    return [
        node
        for node, _ in walk_tree(root)
        if clean_label(str(node.get("name", ""))).casefold() == expected
    ]


def collect_kos(node: dict[str, Any]) -> dict[str, str]:
    """Return KO -> full BRITE entry name for all descendants."""
    result: dict[str, str] = {}

    for child, _ in walk_tree(node):
        name = str(child.get("name", ""))
        match = KO_PATTERN.search(name)

        if match:
            ko_id = match.group(0)
            result[ko_id] = name

    return result


def read_ko_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None or "knum" not in reader.fieldnames:
            raise ValueError(f"{path} does not contain the expected 'knum' column")

        return {
            row["knum"].strip(): {
                key: (value or "").strip() for key, value in row.items()
            }
            for row in reader
            if row.get("knum")
        }


def main() -> None:
    brite_json = Path(snakemake.input.brite_json)
    ko_list = Path(snakemake.input.ko_list)
    profiles_dir = Path(snakemake.params.profiles_dir)

    output_table = Path(snakemake.output.ko_table)
    output_hal = Path(snakemake.output.hal)

    require_thresholds = bool(snakemake.params.require_thresholds)

    with brite_json.open(encoding="utf-8") as handle:
        brite = json.load(handle)

    eukaryote_nodes = find_nodes(brite, "Eukaryotes")

    if len(eukaryote_nodes) != 1:
        raise RuntimeError(
            "Expected exactly one 'Eukaryotes' node in ko03011; "
            f"found {len(eukaryote_nodes)}"
        )

    eukaryote_node = eukaryote_nodes[0]

    subunit_nodes: dict[str, dict[str, Any]] = {}

    for subunit in ("Large subunit", "Small subunit"):
        matches = find_nodes(eukaryote_node, subunit)

        if len(matches) != 1:
            raise RuntimeError(
                f"Expected exactly one '{subunit}' node below "
                f"'Eukaryotes'; found {len(matches)}"
            )

        subunit_nodes[subunit] = matches[0]

    selected: dict[str, dict[str, str]] = {}

    for subunit, node in subunit_nodes.items():
        for ko_id, brite_entry in collect_kos(node).items():
            previous = selected.get(ko_id)

            if previous and previous["subunit"] != subunit:
                raise RuntimeError(f"{ko_id} occurs in more than one subunit")

            selected[ko_id] = {
                "subunit": subunit,
                "brite_entry": brite_entry,
            }

    if not selected:
        raise RuntimeError("No KO identifiers extracted from ko03011")

    ko_metadata = read_ko_metadata(ko_list)

    rows: list[dict[str, str]] = []
    missing_profiles: list[str] = []
    missing_metadata: list[str] = []
    missing_thresholds: list[str] = []

    for ko_id in sorted(selected):
        metadata = ko_metadata.get(ko_id)

        if metadata is None:
            missing_metadata.append(ko_id)
            continue

        profile_path = profiles_dir / f"{ko_id}.hmm"

        if not profile_path.is_file():
            missing_profiles.append(str(profile_path))

        threshold = metadata.get("threshold", "")
        threshold_available = threshold not in {
            "",
            "-",
            "NA",
            "N/A",
        }

        if not threshold_available:
            missing_thresholds.append(ko_id)

        rows.append(
            {
                "ko": ko_id,
                "subunit": selected[ko_id]["subunit"],
                "brite_entry": selected[ko_id]["brite_entry"],
                "threshold": threshold,
                "threshold_available": str(threshold_available).lower(),
                "score_type": metadata.get("score_type", ""),
                "profile_type": metadata.get(
                    "profile_type",
                    "",
                ),
                "definition": metadata.get("definition", ""),
            }
        )

    if missing_metadata:
        raise RuntimeError(
            "KO identifiers absent from local ko_list: " + ", ".join(missing_metadata)
        )

    if missing_profiles:
        raise FileNotFoundError(
            "Missing KOfam profile files:\n" + "\n".join(missing_profiles)
        )

    if require_thresholds and missing_thresholds:
        raise RuntimeError(
            "The following ribosomal KO identifiers do not have "
            "a predefined KOfam threshold: " + ", ".join(missing_thresholds)
        )

    output_table.parent.mkdir(parents=True, exist_ok=True)
    output_hal.parent.mkdir(parents=True, exist_ok=True)

    temporary_table = output_table.with_suffix(".tsv.tmp")

    with temporary_table.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        fieldnames = [
            "ko",
            "subunit",
            "brite_entry",
            "threshold",
            "threshold_available",
            "score_type",
            "profile_type",
            "definition",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)

    temporary_table.replace(output_table)

    hal_lines = [
        "# Cytosolic eukaryotic ribosomal KOfam profiles",
        "# Source: KEGG BRITE ko03011",
        f"# Number of profiles: {len(rows)}",
    ]

    for row in rows:
        profile_path = profiles_dir / f"{row['ko']}.hmm"

        relative_path = os.path.relpath(
            profile_path.resolve(),
            output_hal.parent.resolve(),
        )

        hal_lines.append(relative_path)

    temporary_hal = output_hal.with_suffix(".hal.tmp")
    temporary_hal.write_text(
        "\n".join(hal_lines) + "\n",
        encoding="utf-8",
    )
    temporary_hal.replace(output_hal)


if __name__ == "__main__":
    main()
