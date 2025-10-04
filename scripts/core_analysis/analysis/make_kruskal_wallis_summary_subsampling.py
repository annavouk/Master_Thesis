"""
Biome-level statistics - Kruskal-Wallis

Parse Kruskal-Wallis subsampling results for biomes and compute summary statistics
per Metric.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

import pandas as pd

from config import (
    KRUSKAL_BIOMES_SUBSAMPLES_TSV,  # input
    KRUSKAL_BIOMES_SUMMARY_TSV,  # output
)
from utils import load_data


def make_summary(subs):
    """Build a summary comparing Kruskal-Wallis subsampling dataset."""
    rows = []

    for metric, g in subs.groupby("Metric"):
        sig = g[g["p"] < 0.05]
        n_sig = len(sig)
        perc_sig = (n_sig / len(g)) * 100 if len(g) > 0 else 0

        rows.append(
            {
                "Metric": metric,
                "Significant_reps": f"{n_sig} / {len(g)}",
                "Percent_significant": perc_sig,
                "Median_H_sig": sig["H"].median() if n_sig > 0 else None,
                "Min_H_sig": sig["H"].min() if n_sig > 0 else None,
                "Max_H_sig": sig["H"].max() if n_sig > 0 else None,
                "Median_p_sig": sig["p"].median() if n_sig > 0 else None,
                "Min_p_sig": sig["p"].min() if n_sig > 0 else None,
                "Max_p_sig": sig["p"].max() if n_sig > 0 else None,
            }
        )

    return pd.DataFrame(rows)


def main():
    # Load subsampling results
    subs = load_data(KRUSKAL_BIOMES_SUBSAMPLES_TSV, filetype="tsv")

    summary = make_summary(subs)

    # Save
    summary.to_csv(KRUSKAL_BIOMES_SUMMARY_TSV, sep="\t", index=False)
    print(f"Saved Kruskal summary to {KRUSKAL_BIOMES_SUMMARY_TSV}")
    print(summary)


if __name__ == "__main__":
    main()
