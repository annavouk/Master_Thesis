"""
Quality check script for genome metadata.
Checks how many genomes pass completeness ≥ 90 and contamination ≤ 5,
based on GTDB CheckM estimates. Does not remove any entries.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pandas as pd
from utils import load_data
from config import COMPACT_METADATA


def convert_quality_columns(df, completeness_col='checkm_completeness', contamination_col='checkm_contamination'):
    """Ensure quality columns are numeric."""
    df[completeness_col] = pd.to_numeric(df[completeness_col], errors='coerce')
    df[contamination_col] = pd.to_numeric(df[contamination_col], errors='coerce')
    return df


def quality_report(df, completeness_threshold=90, contamination_threshold=5,
                   completeness_col='checkm_completeness', contamination_col='checkm_contamination'):
    """Print summary of genomes passing quality thresholds."""
    total = len(df)
    passing = df[
        (df[completeness_col] >= completeness_threshold) &
        (df[contamination_col] <= contamination_threshold)
    ]
    n_pass = len(passing)
    percent_pass = (n_pass / total * 100) if total > 0 else 0

    print("Quality Control Report")
    print(f"Total genomes: {total}")
    print(f"Passing genomes: {n_pass} ({percent_pass:.2f}%)")
    print(f"Failing genomes: {total - n_pass}")
    print()
    print("Completeness stats:")
    print(df[completeness_col].describe())
    print()
    print("Contamination stats:")
    print(df[contamination_col].describe())
    print()


def main():
    metadata_path = COMPACT_METADATA

    df = load_data(metadata_path)
    df = convert_quality_columns(df)
    quality_report(df)


if __name__ == "__main__":
    main()
