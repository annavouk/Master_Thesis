"""
Metabolic Interaction Analysis

For every pair of genomes, computes:
1. Absolute feeding matrix (number of seed compounds covered).
2. Normalized feeding matrix (fraction of recipient seeds covered by providers non-seeds).

Generates:
- CSV files of feeding matrices (absolute and normalized)
- Edge lists for network analysis (thresholded)
- Heatmap visualization of top providers/receivers

Usage:
    python3 scripts/analysis/provider_receiver.py --norm-threshold 0.1 --raw-threshold 1 --plot --save-plot --subset 200
"""
import sys
from pathlib import Path
import argparse

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from config import SEEDS_PICKLE, NON_SEEDS_PICKLE, OUTPUT_DIR
from utils import load_data


# ------------------------
# Feeding matrix
# ------------------------
def compute_feeding_matrix(seeds_df, non_seeds_df):
    """
    Return (raw, normalized) feeding matrices.
    Raw = Number of A's seeds that are B's non-seeds
    Normalized = Raw / A's total number of seeds.
    """
    common_cpds = seeds_df.columns.intersection(non_seeds_df.columns)
    seeds_bin = seeds_df[common_cpds].values
    nonseeds_bin = non_seeds_df[common_cpds].values

    abs_mat = np.dot(nonseeds_bin, seeds_bin.T)
    total_seeds = seeds_bin.sum(axis=1)
    total_seeds[total_seeds == 0] = 1
    norm_mat = abs_mat / total_seeds[np.newaxis, :]
    return abs_mat, norm_mat


def feeding_matrix_to_edgelist(feeding_df, threshold):
    """Convert feeding matrix (DataFrame) to edge list for network export."""
    feeding_df = feeding_df.copy()
    feeding_df.index.name = "provider"
    feeding_df.columns.name = "receiver"

    edges = (
        feeding_df
        .stack()
        .reset_index(name="score")
    )
    edges = edges[
        (edges["score"] >= threshold) &
        (edges["provider"] != edges["receiver"])
    ]
    return edges


def plot_feeding_heatmap(feeding_df, top_n=30, save_path=None, show=False):
    """Plot heatmap of top N providers/receivers."""
    top_providers = feeding_df.sum(axis=1).nlargest(top_n).index
    top_receivers = feeding_df.sum(axis=0).nlargest(top_n).index
    sub = feeding_df.loc[top_providers, top_receivers]

    fig, ax = plt.subplots(figsize=(12, 10))
    cax = ax.matshow(sub, cmap="viridis")
    fig.colorbar(cax)

    ax.set_xticks(np.arange(len(top_receivers)))
    ax.set_xticklabels(top_receivers, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(top_providers)))
    ax.set_yticklabels(top_providers)

    ax.set_title("Feeding Matrix Heatmap (Top Providers/Receivers)")
    ax.set_xlabel("Recipient Genome")
    ax.set_ylabel("Provider Genome")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Heatmap saved to {save_path}")
    if show:
        plt.show()
    plt.close(fig)


def print_top_providers_receivers(feeding_df, top_n=10):
    """Print top providers and receivers."""
    top_providers = feeding_df.sum(axis=1).nlargest(top_n)
    top_receivers = feeding_df.sum(axis=0).nlargest(top_n)

    print("\nTop providers (can feed most others):")
    print(top_providers.to_string())
    print("\nTop receivers (most needs met by others):")
    print(top_receivers.to_string())


# ------------------------
# Main
# ------------------------
def main():
    parser = argparse.ArgumentParser(description="Metabolic Interaction Analysis (Feeding Matrix)")
    parser.add_argument("--norm-threshold", type=float, default=0.1,
                        help="Threshold for normalized feeding edge list (default=0.1)")
    parser.add_argument("--raw-threshold", type=int, default=1,
                        help="Threshold for raw feeding edge list (default=1)")
    parser.add_argument("--plot", action="store_true",
                        help="Show heatmap plot interactively")
    parser.add_argument("--save-plot", action="store_true",
                        help="Save heatmap plot as PNG in OUTPUT_DIR")
    parser.add_argument("--subset", type=int, default=None,
                        help="Subset number of genomes for testing (default: all)")
    args = parser.parse_args()

    # Load Pickle dataframes of binary matrices
    seeds_df = load_data(SEEDS_PICKLE, filetype="pickle")
    non_seeds_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")

    if args.subset:
        seeds_df = seeds_df.head(args.subset)
        non_seeds_df = non_seeds_df.head(args.subset)

    # Compute matrices (arrays)
    abs_mat, norm_mat = compute_feeding_matrix(seeds_df, non_seeds_df)

    # Convert to DataFrames
    feeding_raw = pd.DataFrame(abs_mat, index=non_seeds_df.index, columns=seeds_df.index)
    feeding_norm = pd.DataFrame(norm_mat, index=non_seeds_df.index, columns=seeds_df.index)

    # Save matrices
    feeding_raw.to_csv(OUTPUT_DIR / "feeding_matrix_raw.csv")
    feeding_norm.to_csv(OUTPUT_DIR / "feeding_matrix_norm.csv")
    print("Feeding matrices saved: feeding_matrix_raw.csv, feeding_matrix_norm.csv")

    # Export edge lists
    feeding_matrix_to_edgelist(feeding_norm, threshold=args.norm_threshold)\
        .to_csv(OUTPUT_DIR / "feeding_edgelist_norm.csv", index=False)
    feeding_matrix_to_edgelist(feeding_raw, threshold=args.raw_threshold)\
        .to_csv(OUTPUT_DIR / "feeding_edgelist_raw.csv", index=False)
    print("Edge lists exported (raw + normalized)")

    # Heatmap
    if args.plot or args.save_plot:
        plot_feeding_heatmap(
            feeding_norm,
            top_n=30,
            save_path=(OUTPUT_DIR / "feeding_norm_heatmap.png") if args.save_plot else None,
            show=args.plot,
        )

    # Print summary
    print_top_providers_receivers(feeding_norm, top_n=10)


if __name__ == "__main__":
    main()
