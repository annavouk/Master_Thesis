"""
Metabolic Interaction Analysis

Modes:
- In-memory (default): computes full feeding matrices in RAM, exports edge lists, heatmap, summary.
- Blockwise (--blockwise): processes large datasets in blocks, saves matrices + edge lists directly to CSVs.

Usage:
    python3 scripts/analysis/provider_receiver.py --subset 200 --plot --save-plot
    python3 scripts/analysis/provider_receiver.py --blockwise --block-size 1000 --out-prefix run1
"""

import sys
from pathlib import Path
import argparse

sys.path.append(str(Path(__file__).resolve().parents[3]))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from config import SEEDS_PICKLE, NON_SEEDS_PICKLE, OUTPUT_DIR, PLOTS_DIR
from utils import load_data, plot_histogram


# ------------------------
# Helpers for output paths
# ------------------------
def make_output_paths(prefix="feeding", blockwise=False, subset=None):
    """Return dict of output paths with suffixes depending on mode."""
    suffix = []
    if blockwise:
        suffix.append("blockwise")
    if subset:
        suffix.append(f"subset{subset}")
    suf = "_" + "_".join(suffix) if suffix else ""

    return {
        "raw_matrix": OUTPUT_DIR / f"{prefix}_matrix_raw{suf}.csv",
        "norm_matrix": OUTPUT_DIR / f"{prefix}_matrix_norm{suf}.csv",
        "raw_edges": OUTPUT_DIR / f"{prefix}_edgelist_raw{suf}.csv",
        "norm_edges": OUTPUT_DIR / f"{prefix}_edgelist_norm{suf}.csv",
        "heatmap": PLOTS_DIR / f"{prefix}_feeding_heatmap{suf}.png",
        "histogram": PLOTS_DIR / f"{prefix}_edgelist_norm_hist{suf}.png",
    }


# ------------------------
# Core computations
# ------------------------
def compute_feeding_matrix(seeds_df, non_seeds_df):
    """Return (raw, normalized) feeding matrices (in-memory)."""
    common_cpds = seeds_df.columns.intersection(non_seeds_df.columns)
    seeds_bin = seeds_df[common_cpds].values
    nonseeds_bin = non_seeds_df[common_cpds].values

    abs_mat = np.dot(nonseeds_bin, seeds_bin.T)
    total_seeds = seeds_bin.sum(axis=1)
    total_seeds[total_seeds == 0] = 1
    norm_mat = abs_mat / total_seeds[np.newaxis, :]
    return abs_mat, norm_mat


def compute_feeding_matrix_blockwise(seeds_df, non_seeds_df, paths, raw_thr=1, norm_thr=0.1, block_size=1000):
    """Compute feeding matrix block-by-block and save directly to CSVs and edge lists."""
    common_cpds = seeds_df.columns.intersection(non_seeds_df.columns)
    seeds_bin = seeds_df[common_cpds].values.astype(np.uint8)
    nonseeds_bin = non_seeds_df[common_cpds].values.astype(np.uint8)

    total_seeds = seeds_bin.sum(axis=1)
    total_seeds[total_seeds == 0] = 1
    n_providers = nonseeds_bin.shape[0]

    # Init CSVs
    pd.DataFrame(columns=seeds_df.index).to_csv(paths["raw_matrix"], index=False)
    pd.DataFrame(columns=seeds_df.index).to_csv(paths["norm_matrix"], index=False)
    pd.DataFrame(columns=["provider", "receiver", "score"]).to_csv(paths["raw_edges"], index=False)
    pd.DataFrame(columns=["provider", "receiver", "score"]).to_csv(paths["norm_edges"], index=False)

    # Process in chunks
    for start in range(0, n_providers, block_size):
        end = min(start + block_size, n_providers)
        block = nonseeds_bin[start:end, :]

        abs_block = np.dot(block, seeds_bin.T)
        norm_block = abs_block / total_seeds[np.newaxis, :]

        providers = non_seeds_df.index[start:end]
        receivers = seeds_df.index

        # Save matrices
        raw_df = pd.DataFrame(abs_block, index=providers, columns=receivers)
        norm_df = pd.DataFrame(norm_block, index=providers, columns=receivers)
        raw_df.to_csv(paths["raw_matrix"], mode="a", header=False)
        norm_df.to_csv(paths["norm_matrix"], mode="a", header=False)

        # Save edges
        raw_edges = (
            raw_df.stack().reset_index(name="score")
            .query("score >= @raw_thr and provider != receiver")
        )
        norm_edges = (
            norm_df.stack().reset_index(name="score")
            .query("score >= @norm_thr and provider != receiver")
        )
        raw_edges.to_csv(paths["raw_edges"], mode="a", header=False, index=False)
        norm_edges.to_csv(paths["norm_edges"], mode="a", header=False, index=False)

        print(f"Processed block {start}:{end} / {n_providers}")

    print(f"\nBlockwise outputs saved:")
    for k, v in paths.items():
        if "heatmap" not in k and "histogram" not in k:
            print(f" - {v}")


# ------------------------
# Visualization helpers
# ------------------------
def feeding_matrix_to_edgelist(feeding_df, threshold):
    """Convert feeding matrix (DataFrame) to edge list for network export."""
    feeding_df = feeding_df.copy()
    feeding_df.index.name = "provider"
    feeding_df.columns.name = "receiver"

    edges = feeding_df.stack().reset_index(name="score")
    edges = edges[(edges["score"] >= threshold) & (edges["provider"] != edges["receiver"])]
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


def plot_top_providers_receivers_from_edges(edges_csv, top_n=20, plots_dir=PLOTS_DIR, prefix="feeding"):
    """Generate barplots of top providers and receivers from an edge list CSV."""
    edges = pd.read_csv(edges_csv)

    provider_strength = edges.groupby("provider")["score"].sum().nlargest(top_n)
    receiver_strength = edges.groupby("receiver")["score"].sum().nlargest(top_n)

    # Providers
    plt.figure(figsize=(10, 6))
    sns.barplot(x=provider_strength.values, y=provider_strength.index, color="steelblue")
    plt.xlabel("Total normalized feeding score")
    plt.ylabel("Provider genome")
    plt.title(f"Top {top_n} Providers")
    plt.tight_layout()
    out_path = plots_dir / f"{prefix}_top{top_n}_providers.png"
    plt.savefig(out_path, dpi=300)
    plt.close()

    # Receivers
    plt.figure(figsize=(10, 6))
    sns.barplot(x=receiver_strength.values, y=receiver_strength.index, color="darkorange")
    plt.xlabel("Total normalized feeding score")
    plt.ylabel("Receiver genome")
    plt.title(f"Top {top_n} Receivers")
    plt.tight_layout()
    out_path = plots_dir / f"{prefix}_top{top_n}_receivers.png"
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"Top {top_n} providers/receivers barplots saved to {plots_dir}")


def print_top_providers_receivers(feeding_df, top_n=10):
    """Print top providers and receivers."""
    top_providers = feeding_df.sum(axis=1).nlargest(top_n)
    top_receivers = feeding_df.sum(axis=0).nlargest(top_n)
    print("\nTop providers:\n", top_providers.to_string())
    print("\nTop receivers:\n", top_receivers.to_string())


# ------------------------
# Run modes
# ------------------------
def run_in_memory(seeds_df, non_seeds_df, args, paths):
    abs_mat, norm_mat = compute_feeding_matrix(seeds_df, non_seeds_df)
    feeding_raw = pd.DataFrame(abs_mat, index=non_seeds_df.index, columns=seeds_df.index)
    feeding_norm = pd.DataFrame(norm_mat, index=non_seeds_df.index, columns=seeds_df.index)

    feeding_raw.to_csv(paths["raw_matrix"])
    feeding_norm.to_csv(paths["norm_matrix"])
    print("Feeding matrices saved (in-memory mode).")

    feeding_matrix_to_edgelist(feeding_norm, threshold=args.norm_threshold)\
        .to_csv(paths["norm_edges"], index=False)
    feeding_matrix_to_edgelist(feeding_raw, threshold=args.raw_threshold)\
        .to_csv(paths["raw_edges"], index=False)
    print("Edge lists exported.")

    if args.plot or args.save_plot:
        plot_feeding_heatmap(
            feeding_norm,
            top_n=30,
            save_path=paths["heatmap"] if args.save_plot else None,
            show=args.plot,
        )

    print_top_providers_receivers(feeding_norm, top_n=10)


def run_blockwise(seeds_df, non_seeds_df, args, paths):
    compute_feeding_matrix_blockwise(
        seeds_df, non_seeds_df,
        paths=paths,
        raw_thr=args.raw_threshold,
        norm_thr=args.norm_threshold,
        block_size=args.block_size
    )
    print("Heatmap not available in blockwise mode (matrix too large).")


# ------------------------
# Main
# ------------------------
def main():
    parser = argparse.ArgumentParser(description="Metabolic Interaction Analysis (Feeding Matrix)")
    parser.add_argument("--norm-threshold", type=float, default=0.1,
                        help="Threshold for normalized feeding edge list (default=0.1)")
    parser.add_argument("--raw-threshold", type=int, default=1,
                        help="Threshold for raw feeding edge list (default=1)")
    parser.add_argument("--plot", action="store_true", help="Show heatmap plot interactively")
    parser.add_argument("--save-plot", action="store_true", help="Save heatmap plot as PNG")
    parser.add_argument("--subset", type=int, default=None, help="Subset of genomes for testing")
    parser.add_argument("--blockwise", action="store_true", help="Use blockwise computation")
    parser.add_argument("--block-size", type=int, default=1000, help="Block size for blockwise mode")
    parser.add_argument("--out-prefix", type=str, default="feeding", help="Prefix for output files")
    parser.add_argument("--histogram", action="store_true", help="Make histogram of normalized edge scores")
    parser.add_argument("--barplots", action="store_true", 
                    help="Generate barplots of top providers/receivers from normalized edge list")
    args = parser.parse_args()

    # Load Pickle dataframes
    seeds_df = load_data(SEEDS_PICKLE, filetype="pickle")
    non_seeds_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")
    if args.subset:
        seeds_df = seeds_df.head(args.subset)
        non_seeds_df = non_seeds_df.head(args.subset)

    paths = make_output_paths(prefix=args.out_prefix, blockwise=args.blockwise, subset=args.subset)

    if args.blockwise:
        run_blockwise(seeds_df, non_seeds_df, args, paths)
    else:
        run_in_memory(seeds_df, non_seeds_df, args, paths)

    # Extra: Histogram from normalized edge list
    if args.histogram:
        edges = pd.read_csv(paths["norm_edges"])
        
        plot_histogram(
            edges,
            column="score",
            bins=50,
            color="teal",
            title="Distribution of normalized feeding scores",
            xlabel="Normalized feeding score",
            save_plot=True,
            plot_path=paths["histogram"],
            log_scale=False,
            method="2std"
        )

    # Extra: Barplots from normalized edge list
    if args.barplots:
        plot_top_providers_receivers_from_edges(
            edges_csv=paths["norm_edges"],
            top_n=20,
            plots_dir=PLOTS_DIR,
            prefix=args.out_prefix
        )


if __name__ == "__main__":
    main()
