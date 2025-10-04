"""
Biome-level statistics - Dunn's post-hoc

Parse Dunn's post-hoc subsampling results for biomes and compute summary statistics
per Metric/Pair combination.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

import pandas as pd
from io import StringIO

from config import (
    DUNN_BIOMES_SUBSAMPLES_TSV,    # input
    DUNN_BIOMES_SUMMARY_TSV,    # output
    )


def make_dunn_summary(subs_path):
    pairs = []

    with open(subs_path) as f:
        rep = None
        metric = None
        block = []

        for line in f:
            if line.startswith("# Rep"):
                if block:  # process previous block
                    df = pd.read_csv(StringIO("".join(block)), sep="\t", index_col=0)
                    melted = df.stack().reset_index()
                    melted.columns = ["Group1", "Group2", "p_value"]
                    melted = melted[melted["Group1"] != melted["Group2"]]
                    melted["Pair"] = melted.apply(
                        lambda r: " vs ".join(sorted([r["Group1"], r["Group2"]])), axis=1
                    )
                    melted["Metric"] = metric
                    melted["Rep"] = rep
                    pairs.append(melted[["Rep", "Metric", "Pair", "p_value"]])
                    block = []

                parts = line.strip().split()
                rep = int(parts[2])
                metric = parts[-1]

            elif line.strip():
                block.append(line)

        # process last block
        if block:
            df = pd.read_csv(StringIO("".join(block)), sep="\t", index_col=0)
            melted = df.stack().reset_index()
            melted.columns = ["Group1", "Group2", "p_value"]
            melted = melted[melted["Group1"] != melted["Group2"]]
            melted["Pair"] = melted.apply(
                lambda r: " vs ".join(sorted([r["Group1"], r["Group2"]])), axis=1
            )
            melted["Metric"] = metric
            melted["Rep"] = rep
            pairs.append(melted[["Rep", "Metric", "Pair", "p_value"]])

    if not pairs:
        raise ValueError("No blocks found! Check input file format.")

    dunn_subs_long = pd.concat(pairs, ignore_index=True)

    # Summary
    summary = (
        dunn_subs_long
        .groupby(["Metric", "Pair"])
        .agg(
            Significant_reps=("p_value", lambda x: f"{(x<0.05).sum()} / {len(x)}"),
            Percent_significant=("p_value", lambda x: 100*(x<0.05).mean()),
            Median_p_sig=("p_value", lambda x: x[x<0.05].median() if any(x<0.05) else None),
            Min_p_sig=("p_value", lambda x: x[x<0.05].min() if any(x<0.05) else None),
            Max_p_sig=("p_value", lambda x: x[x<0.05].max() if any(x<0.05) else None),
        )
        .reset_index()
    )

    return summary


def main():
    # Load data
    subs_path = Path(DUNN_BIOMES_SUBSAMPLES_TSV)

    # Compute summary per Metric/Pair
    subs_summary = make_dunn_summary(subs_path)

    # Save
    subs_summary.to_csv(DUNN_BIOMES_SUMMARY_TSV, sep="\t", index=False)
    print(f"Saved Dunn's summary (subs) to {DUNN_BIOMES_SUMMARY_TSV}")
    print(subs_summary.head())


if __name__ == "__main__":
    main()
