#!/usr/bin/env python


"""
Program reads table from tRNAscan-SE output and count previously complemented codons. Output is in format needed in gtAI (R).

Note: input Anti Codon is in DNA (as seq in input tRNAscan-SE)


!IMPORTANT!

in the output's column 'Anti codon' there are sequences of tRNA anticodons, f.e. AAT in DNA 5'-> 3'. They are not codons from coding sequence.

tRNAscan-SE reports tRNA-Ile with AAT anticodon, so biological-real-life sequence in tRNA is RNA 5'-AAU-3'. Proper codon in CDS is reverse complemented: DNA 5'-ATT-3' and that

Input:
    - tRNAscanSE.out file

Output:
    - .tsv

Usage:


Example:

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
        prog="count_codons_from_trnascan.py",
        description="Count codons occurences from tRNAscan-SE output file converted to .csv file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        conflict_handler="error",
        add_help=True,
    )

    p.add_argument(
        "-I",
        "--input",
        help="Path to input file",
        required=True,
        type=Path,
    )

    p.add_argument(
        "-O",
        "--output",
        help="Path to ouptut file",
        required=True,
        type=Path,
    )


# ------------------------
# --- Helper functions ---
# ------------------------

# ---------------------
# --- Main function ---
# ---------------------
