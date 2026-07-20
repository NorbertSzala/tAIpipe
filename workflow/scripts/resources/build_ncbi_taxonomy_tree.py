#!/usr/bin/env python3
"""Build a low-cost NCBI taxonomy topology from genome-accession metadata.

Input is the TSV produced by:

    datasets summary genome accession --inputfile accessions.txt --as-json-lines |
      dataformat tsv genome \
        --fields accession,organism-name,organism-tax-id \
        > accession_taxonomy.tsv

The result is a taxonomy-derived topology, not a sequence-inferred species tree.
When several accessions share one TaxID, they are added as zero-length sister
leaves under that taxon because their within-species relationship is unknown.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def normalize_header(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower())


def choose_column(fieldnames: list[str], candidates: set[str], label: str) -> str:
    normalized = {normalize_header(name): name for name in fieldnames}
    for candidate in candidates:
        if candidate in normalized:
            return normalized[candidate]
    raise ValueError(f"Could not identify the {label} column in: {fieldnames}")


def safe_label(value: str) -> str:
    value = value.strip()
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"[():;,\[\]'\"]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--leaf-label",
        choices=("accession", "organism_accession"),
        default="accession",
    )
    parser.add_argument(
        "--update-taxonomy",
        action="store_true",
        help="Refresh the local ETE NCBI taxonomy database before building the tree.",
    )
    args = parser.parse_args()

    try:
        from ete3 import NCBITaxa
    except ImportError as exc:
        raise SystemExit(
            "ETE3 is required. Install it, for example, with: "
            "mamba install -c conda-forge -c bioconda ete3"
        ) from exc

    with args.metadata.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Metadata TSV has no header")

        accession_col = choose_column(
            reader.fieldnames,
            {"accession", "assemblyaccession"},
            "accession",
        )
        organism_col = choose_column(
            reader.fieldnames,
            {"organismname", "organism"},
            "organism name",
        )
        taxid_col = choose_column(
            reader.fieldnames,
            {"organismtaxid", "organismtaxonomicid", "taxid"},
            "TaxID",
        )

        records: list[dict[str, str | int]] = []
        for row in reader:
            accession = row.get(accession_col, "").strip()
            organism = row.get(organism_col, "").strip()
            taxid_raw = row.get(taxid_col, "").strip()
            if not accession or not taxid_raw:
                continue
            try:
                taxid = int(taxid_raw)
            except ValueError as exc:
                raise ValueError(f"Invalid TaxID {taxid_raw!r} for {accession}") from exc
            records.append(
                {"accession": accession, "organism": organism, "taxid": taxid}
            )

    if not records:
        raise ValueError("No usable accession/TaxID records were found")

    by_taxid: dict[int, list[dict[str, str | int]]] = defaultdict(list)
    for record in records:
        by_taxid[int(record["taxid"])].append(record)

    ncbi = NCBITaxa()
    if args.update_taxonomy:
        ncbi.update_taxonomy_database()

    taxids = sorted(by_taxid)
    tree = ncbi.get_topology(taxids, intermediate_nodes=True)

    # Unit edges encode taxonomic/topological steps only; they are not
    # divergence times. Set them before adding zero-length replicate leaves.
    for node in tree.traverse():
        if not node.is_root():
            node.dist = 1.0

    # Add the requested accession labels before renaming taxonomic nodes.
    for taxid, members in by_taxid.items():
        matches = tree.search_nodes(name=str(taxid))
        if not matches:
            raise RuntimeError(f"TaxID {taxid} is absent from the extracted topology")
        node = matches[0]

        labels = []
        for member in members:
            accession = str(member["accession"])
            organism = str(member["organism"])
            if args.leaf_label == "organism_accession":
                labels.append(safe_label(f"{organism}_{accession}"))
            else:
                labels.append(safe_label(accession))

        if len(labels) == 1 and node.is_leaf():
            node.name = labels[0]
        else:
            # The taxonomy does not resolve relationships among assemblies that
            # share a TaxID. Represent them as an explicit zero-length polytomy.
            for label in labels:
                node.add_child(name=label, dist=0.0)

    numeric_node_ids = [
        int(node.name)
        for node in tree.traverse()
        if node.name and node.name.isdigit()
    ]
    taxid_to_name = ncbi.get_taxid_translator(numeric_node_ids)
    for node in tree.traverse():
        if node.name and node.name.isdigit():
            node.name = safe_label(taxid_to_name.get(int(node.name), node.name))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    tree.write(outfile=str(args.output), format=1)
    print(f"Wrote taxonomy topology for {len(records)} accessions to {args.output}")
    duplicated = sum(len(v) > 1 for v in by_taxid.values())
    if duplicated:
        print(
            f"Warning: {duplicated} TaxIDs contained multiple accessions; "
            "those accessions were represented as unresolved sister tips."
        )


if __name__ == "__main__":
    main()