#!/usr/bin/env python3


'''
## Description

`parse_go_obo_to_tsv.py` creates a GO term dictionary TSV from a Gene Ontology OBO file.

If a valid output TSV already exists, the script exits without changes. Otherwise, it parses the provided OBO file (or downloads it if missing and `--obo-url` is supplied) and generates a table containing GO IDs, names, and namespaces.

## Input

Required:
- `--input-obo` -  path to the GO `.obo` file.
- `--output-tsv` – output GO dictionary TSV.

Optional:
- `--obo-url` – URL used to download the OBO file when missing.

## Output

TSV file with columns:

- `go_terms`
- `go_name`
- `go_namespace`
'''

# ----------------------
# --- Import modules ---
# ----------------------

from pathlib import Path
import argparse
import urllib.request
import sys


# -----------------------
# --- Constant values ---
# -----------------------
REQUIRED_COLUMNS = ["go_terms", "go_name", "go_namespace"]


# -----------------
# --- Arguments ---
# -----------------
def arguments():
    p = argparse.ArgumentParser(
        prog = "parse_go_obo_to_tsv.py",
        description="Prepare GO term dictionary TSV from GO OBO file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        conflict_handler="error",
        add_help=True,
    )
    
    p.add_argument("-I", "--input-obo", required=True, type = Path)
    p.add_argument("-O", "--output-tsv", required=True, type=Path)
    p.add_argument("-U", "--obo-url", type = str, default="")
    return p.parse_args()


# ------------------------
# --- Helper functions ---
# ------------------------
def tsv_is_valid(path: Path) -> bool:
    """Checks whether an existing tsv has the required columns"""
    if not path.exists() or path.stat().st_size == 0:
        return False

    with path.open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")

    return all(col in header for col in REQUIRED_COLUMNS)


def download_obo(url: str, output_path: Path) -> None:
    """Download GO OBO file to the requested local path"""
    if not url:
        raise ValueError("OBO file is missing and no --obo-url was provided.")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Downloading GO OBO from: {url}", file=sys.stderr)
    print(f"Saving to: {output_path}", file=sys.stderr)

    urllib.request.urlretrieve(url, output_path)


def parse_go_obo(input_obo: Path, output_tsv: Path) -> None:
    rows = []
    current = {}

    with input_obo.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if line == "[Term]":
                if (
                    current.get("id")
                    and current.get("name")
                    and current.get("namespace")
                    and not current.get("is_obsolete", False)
                ):
                    rows.append(current)

                current = {}
                continue

            if not line or line.startswith("!"):
                continue

            if line.startswith("id: "):
                current["id"] = line.replace("id: ", "", 1)

            elif line.startswith("name: "):
                current["name"] = line.replace("name: ", "", 1)

            elif line.startswith("namespace: "):
                current["namespace"] = line.replace("namespace: ", "", 1)

            elif line.startswith("is_obsolete: true"):
                current["is_obsolete"] = True

    if (
        current.get("id")
        and current.get("name")
        and current.get("namespace")
        and not current.get("is_obsolete", False)
    ):
        rows.append(current)

    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    with output_tsv.open("w", encoding="utf-8") as out:
        out.write("go_terms\tgo_name\tgo_namespace\n")
        for row in rows:
            out.write(f"{row['id']}\t{row['name']}\t{row['namespace']}\n")

    print(f"Written {len(rows)} GO terms to: {output_tsv}", file=sys.stderr)


# ---------------------
# --- Main function ---
# ---------------------

def main() -> None:
    "Prepare go dictionary tsv from existing tsv, local obo or downloaded OBO"
    args = arguments()
    input_obo = args.input_obo
    output_tsv = args.output_tsv


    if tsv_is_valid(output_tsv):
        print(f"Existing valid GO dictionary TSV found: {output_tsv}", file=sys.stderr)
        return

    if not input_obo.exists():
        print(f"OBO file not found: {input_obo}", file = sys.stderr)
        download_obo(args.obo_url, input_obo)

    if not input_obo.exists() or input_obo.stat().st_size == 0:
        raise FileNotFoundError(f"OBO file not found or empty: {input_obo}")

    parse_go_obo(input_obo, output_tsv)

    if not tsv_is_valid(output_tsv):
        raise RuntimeError(f"Generated TSV is invalid: {output_tsv}")


if __name__ == "__main__":
    main()