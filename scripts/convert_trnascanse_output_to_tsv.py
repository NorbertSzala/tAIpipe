#!/usr/bin/env python

"""
Description:

    tRNAscan-SE 2.0 output is very awful:

    Sequence   		tRNA    	Bounds  	tRNA	Anti	Intron Bounds	Inf
    Name       	tRNA #	Begin   	End     	Type	Codon	Begin	End	Score	Note
    --------   	------	-----   	------  	----	-----	-----	----	------	------
    LR898468.1 	1	656159  	656232  	Ile	AAT	0	0	78.2
    LR898468.1 	2	701452  	701523  	Arg	TCG	0	0	69.1
    LR898468.1 	3	837869  	837965  	Leu	AAG	837907	837921	40.7
    LR898468.1 	4	1085856 	1085927 	Trp	CCA	0	0	65.5
    LR898468.1 	5	1121225 	1121306 	Asp	GTC	1121262	1121271	54.8


    The aim of this script is to convert and validate the given file into real .tsv format to make future analysis much easier.


Input:
    - tRNAscanSE.out file

Output:
    - clean .tsv

    output columns:
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
Usage:


Example usage:

    ./scripts/convert_trnascanse_output_to_tsv.py -I results/per_genome/Spombe/trnascan/Spombe_trnascan.out -O ./try.tsv --remove-NNN
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
        prog="convert_trnascanse_output_to_tsv.py",
        description="Convert tRNAscan-SE output into clean TSV format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        conflict_handler="error",
        add_help=True,
    )

    p.add_argument(
        "-I",
        "--input",
        help="Path to input tRNAscan-SE .out file",
        required=True,
        type=Path,
    )

    p.add_argument(
        "-O",
        "--output",
        help="Path to output .tsv file",
        required=True,
        type=Path,
    )

    p.add_argument(
        "--keep-pseudo",
        help="Keep tRNAs marked as pseudo",
        action="store_true",
    )

    p.add_argument(
        "--remove-NNN",
        help="Remove tRNAs with unrecognized anticodon written as NNN in tRNAscan-SE output",
        action="store_true",
    )

    return p.parse_args()


# ------------------------
# --- Helper functions ---
# ------------------------


def validate_input_file(input_path: Path) -> None:
    """Check whether input file exists and is not empty."""
    if not input_path.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_path}")

    if not input_path.is_file():
        raise ValueError(f"Input path is not a file: {input_path}")

    if input_path.stat().st_size == 0:
        raise ValueError(f"Input file is empty: {input_path}")

    if input_path.suffix != ".out":
        raise Warning(
            f"Input file has improper extenstion: {input_path.suffix}. It should be '.out'"
        )


def is_data_line(line: str) -> bool:
    """
    Return True only for real data lines.

    Skips:
    - empty lines
    - header lines
    - separator lines
    - comment-like lines
    """
    line = line.strip()

    if not line:
        return False

    if line.startswith("-"):
        return False

    if line.startswith("Sequence"):
        return False

    if line.startswith("Name"):
        return False

    if line.startswith("Score"):
        return False

    return True


def parse_trnascanse_output(input_path: Path) -> pd.DataFrame:
    """
    Parse tRNAscan-SE output into pandas DataFrame.

    Expected data columns:
    seq_name, trna_number, begin, end, trna_type, anticodon,
    intron_begin, intron_end, score, optional note
    """
    records = []

    with input_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not is_data_line(line):
                continue

            fields = line.strip().split()

            if len(fields) not in [9, 10]:
                raise ValueError(
                    f"Unexpected number of columns in line {line_number}: "
                    f"expected 9 or 10, got {len(fields)}\n"
                    f"Line content: {line.rstrip()}"
                )

            if len(fields) == 9:
                fields.append("")

            (
                seq_name,
                trna_number,
                begin,
                end,
                trna_type,
                anticodon,
                intron_begin,
                intron_end,
                score,
                note,
            ) = fields

            records.append(
                {
                    "seq_name": seq_name,
                    "trna_number": int(trna_number),
                    "begin": int(begin),
                    "end": int(end),
                    "trna_type": trna_type,
                    "anticodon": anticodon,
                    "intron_begin": int(intron_begin),
                    "intron_end": int(intron_end),
                    "score": float(score),
                    "note": note,
                }
            )

    if not records:
        raise ValueError(f"No valid tRNAscan-SE records found in: {input_path}")

    df = pd.DataFrame(records)

    return df


def add_derived_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add useful columns:
    - strand: '+' if begin <= end, otherwise '-'
    - start: smaller coordinate
    - stop: larger coordinate - proper for 5' -> 3' format
    - has_intron: True if intron coordinates are not 0
    - is_pseudo: True if note contains 'pseudo'
    """
    df = df.copy()

    df["strand"] = df.apply(
        lambda row: "+" if row["begin"] <= row["end"] else "-",
        axis=1,
    )

    df["start"] = df[["begin", "end"]].min(axis=1)
    df["stop"] = df[["begin", "end"]].max(axis=1)

    df["has_intron"] = (df["intron_begin"] != 0) & (df["intron_end"] != 0)

    df["is_pseudo"] = (
        df["note"].fillna("").str.lower().str.contains("pseudo", regex=False)
    )

    return df


def validate_dataframe(df: pd.DataFrame) -> None:
    """Basic sanity checks."""
    required_columns = [
        "seq_name",
        "trna_number",
        "begin",
        "end",
        "trna_type",
        "anticodon",
        "intron_begin",
        "intron_end",
        "score",
        "note",
        "strand",
        "start",
        "stop",
        "has_intron",
        "is_pseudo",
    ]

    missing = set(required_columns) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    invalid_anticodons = df[~df["anticodon"].str.fullmatch(r"[ACGTN]{3}", na=False)]

    if not invalid_anticodons.empty:
        raise ValueError(
            "Invalid anticodon values found:\n"
            f"{invalid_anticodons[['seq_name', 'trna_number', 'anticodon']]}"
        )

    if (df["score"] < 0).any():
        raise ValueError("Negative Infernal scores found; this is unexpected.")


def write_output(df: pd.DataFrame, output_path: Path) -> None:
    """Write DataFrame to TSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(
        output_path,
        sep="\t",
        index=False,
    )


# ---------------------
# --- Main function ---
# ---------------------
def main():
    args = arguments()

    validate_input_file(args.input)

    df = parse_trnascanse_output(args.input)
    df = add_derived_columns(df)
    validate_dataframe(df)

    if not args.keep_pseudo:
        df = df[~df["is_pseudo"]].copy()

    if args.remove_NNN:
        df = df[df["anticodon"] != "NNN"].copy()

    write_output(df, args.output)

    print(f"Saved clean TSV: {args.output}")
    print(f"Records written: {len(df)}")


if __name__ == "__main__":
    main()
