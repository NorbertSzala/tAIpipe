#!/usr/bin/env python3

"""
# Creates a compact subset of the Gene Ontology OBO file for tests. It preserves the GO records needed by the test dataset so that ontology parsing and enrichment steps can be exercised without downloading or processing the complete ontology.

Create a small Gene Ontology OBO subset for tests.

The script:
1. Reads a TSV file containing a semicolon-separated GO term column.
2. Extracts all unique GO IDs used in that file.
3. Reads a full GO OBO file.
4. Keeps only [Term] records matching the extracted GO IDs.
5. Writes a small test OBO file suitable for committing to the repository.
"""

from pathlib import Path
import argparse
import csv
import sys


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog="create_test_go_obo_subset.py",
        description="Extract GO terms from a TSV file and create a matching small OBO subset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-t",
        "--input-tsv",
        required=True,
        type=Path,
        help="Input TSV file containing a GO terms column, e.g. samples_test.tsv.",
    )

    parser.add_argument(
        "-g",
        "--go-column",
        default="go_terms",
        help="Name of the column containing semicolon-separated GO IDs.",
    )

    parser.add_argument(
        "-o",
        "--input-obo",
        required=True,
        type=Path,
        help="Full Gene Ontology OBO file.",
    )

    parser.add_argument(
        "-O",
        "--output-obo",
        required=True,
        type=Path,
        help="Output small OBO subset file.",
    )

    parser.add_argument(
        "--output-terms",
        type=Path,
        default=None,
        help="Optional output text file with extracted GO IDs, one per line.",
    )

    return parser.parse_args()


def extract_go_terms_from_tsv(input_tsv: Path, go_column: str) -> set[str]:
    """Extract unique GO IDs from a semicolon-separated GO column."""
    if not input_tsv.exists():
        raise FileNotFoundError(f"Input TSV file does not exist: {input_tsv}")

    go_terms = set()

    with input_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        if reader.fieldnames is None:
            raise ValueError(f"Input TSV has no header: {input_tsv}")

        if go_column not in reader.fieldnames:
            raise ValueError(
                f"Column '{go_column}' not found in {input_tsv}. "
                f"Available columns: {', '.join(reader.fieldnames)}"
            )

        for row in reader:
            raw_value = row.get(go_column, "")

            if raw_value is None:
                continue

            raw_value = raw_value.strip()

            if not raw_value:
                continue

            for term in raw_value.split(";"):
                term = term.strip()

                if term:
                    go_terms.add(term)

    return go_terms


def split_obo_into_header_and_terms(
    input_obo: Path,
) -> tuple[list[str], list[list[str]]]:
    """Split an OBO file into header lines and individual [Term] blocks."""
    if not input_obo.exists():
        raise FileNotFoundError(f"Input OBO file does not exist: {input_obo}")

    header_lines = []
    term_blocks = []
    current_block = None

    with input_obo.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()

            if stripped == "[Term]":
                if current_block is not None:
                    term_blocks.append(current_block)

                current_block = [line]
                continue

            if current_block is None:
                header_lines.append(line)
            else:
                current_block.append(line)

    if current_block is not None:
        term_blocks.append(current_block)

    return header_lines, term_blocks


def get_term_id(term_block: list[str]) -> str | None:
    """Return GO ID from a single [Term] block."""
    for line in term_block:
        if line.startswith("id: "):
            return line.replace("id: ", "", 1).strip()

    return None


def is_obsolete(term_block: list[str]) -> bool:
    """Check whether a [Term] block is marked as obsolete."""
    return any(line.strip() == "is_obsolete: true" for line in term_block)


def write_go_terms(go_terms: set[str], output_terms: Path) -> None:
    """Write extracted GO IDs to a plain text file."""
    output_terms.parent.mkdir(parents=True, exist_ok=True)

    with output_terms.open("w", encoding="utf-8") as out:
        for term in sorted(go_terms):
            out.write(term + "\n")


def write_obo_subset(
    header_lines: list[str],
    selected_blocks: list[list[str]],
    output_obo: Path,
) -> None:
    """Write selected OBO [Term] blocks to a new OBO file."""
    output_obo.parent.mkdir(parents=True, exist_ok=True)

    with output_obo.open("w", encoding="utf-8") as out:
        if header_lines:
            out.writelines(header_lines)

            if header_lines[-1].strip() != "":
                out.write("\n")
        else:
            out.write("format-version: 1.2\n")
            out.write("data-version: test-go-subset\n\n")

        for block in selected_blocks:
            if block and block[-1].strip() != "":
                block = block + ["\n"]

            out.writelines(block)


def main() -> None:
    args = parse_args()

    go_terms = extract_go_terms_from_tsv(
        input_tsv=args.input_tsv,
        go_column=args.go_column,
    )

    if not go_terms:
        raise ValueError(
            f"No GO terms were found in column '{args.go_column}' "
            f"of file: {args.input_tsv}"
        )

    header_lines, term_blocks = split_obo_into_header_and_terms(args.input_obo)

    selected_blocks = []
    found_terms = set()

    for block in term_blocks:
        term_id = get_term_id(block)

        if term_id is None:
            continue

        if term_id in go_terms and not is_obsolete(block):
            selected_blocks.append(block)
            found_terms.add(term_id)

    missing_terms = sorted(go_terms - found_terms)

    write_obo_subset(
        header_lines=header_lines,
        selected_blocks=selected_blocks,
        output_obo=args.output_obo,
    )

    if args.output_terms is not None:
        write_go_terms(go_terms, args.output_terms)

    print(f"Input TSV: {args.input_tsv}", file=sys.stderr)
    print(f"Input OBO: {args.input_obo}", file=sys.stderr)
    print(f"Extracted GO terms from TSV: {len(go_terms)}", file=sys.stderr)
    print(f"GO terms found in OBO: {len(found_terms)}", file=sys.stderr)
    print(f"GO terms missing in OBO: {len(missing_terms)}", file=sys.stderr)
    print(f"Written OBO subset: {args.output_obo}", file=sys.stderr)

    if args.output_terms is not None:
        print(f"Written GO term list: {args.output_terms}", file=sys.stderr)

    if missing_terms:
        print(
            "\nWarning: some GO terms from the TSV were not found in the OBO file:",
            file=sys.stderr,
        )
        for term in missing_terms[:50]:
            print(f"  {term}", file=sys.stderr)

        if len(missing_terms) > 50:
            print(f"  ... and {len(missing_terms) - 50} more", file=sys.stderr)


if __name__ == "__main__":
    main()
