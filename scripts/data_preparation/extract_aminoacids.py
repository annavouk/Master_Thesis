"""
Extract the 20 common amino acids from the compound summary file.

- Loads compound_summary.tsv
- Filters compounds based on KEGG IDs of standard amino acids as described in KEGG Compounds with Biological Roles br080001
(Peptides>Peptides>Amino acids>Common amino acids)
- Prints and saves the subset to AMINOACIDS_TSV
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd

from config import (
    COMPOUND_SUMMARY_TSV,  # input
    AMINOACIDS_TSV,  # output
)
from utils import load_data


# ------------------------
# Extract common aminoacids
# ------------------------
def get_aa(df):
    """Filter compound summary for the 20 standard amino acids (KEGG IDs)."""
    common_aa_kegg_ids = [
        "C00037",
        "C00041",
        "C00183",
        "C00123",
        "C00407",
        "C00049",
        "C00152",
        "C00025",
        "C00064",
        "C00065",
        "C00188",
        "C00073",
        "C00097",
        "C00047",
        "C00062",
        "C00135",
        "C00148",
        "C00079",
        "C00082",
        "C00078",
    ]

    common_aa_df = df[df["KEGG_ID"].isin(common_aa_kegg_ids)]

    print(f"Found {len(common_aa_df)} compounds matching common amino acids.")
    print(common_aa_df[["ModelSEED_ID", "compound_name", "KEGG_ID"]])

    return common_aa_df


# ------------------------
# Main
# ------------------------
def main():
    # Load data
    df = load_data(COMPOUND_SUMMARY_TSV, filetype="tsv")

    # Extract aminoacids
    common_aa_df = get_aa(df)

    # Save output
    common_aa_df.to_csv(AMINOACIDS_TSV, sep="\t", index=False)
    print(f"Aminoacids saved to {AMINOACIDS_TSV}")


if __name__ == "__main__":
    main()
