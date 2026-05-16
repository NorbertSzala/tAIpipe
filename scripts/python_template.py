#!/usr/bin/env python

"""

Description:

Input:

Output:

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
        prog="",
        description="",
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

    return p.parse_args()


# ------------------------
# --- Helper functions ---
# ------------------------

# ---------------------
# --- Main function ---
# ---------------------
