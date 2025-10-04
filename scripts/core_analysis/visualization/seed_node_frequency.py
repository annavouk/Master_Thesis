"""
Seed node frequency

Analyze the frequency of seed nodes (essential metabolites) across genomes.
    -Histogram: Distribution of seed frequencies across compounds.
    -Barplot: Top-N most frequent seed compounds.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from config import (
    SEEDS_PICKLE,  # input
    NON_SEEDS_PICKLE,  # input
    COMPOUND_SUMMARY_TSV,  # metadata
    SEED_NODE_FREQUENCY_PLOTS_DIR,  # plots dir
    OUTPUT_DIR,  # output dir
)
from utils import load_data


# ------------------------
# Compute seed node frequency
# ------------------------
def compute_seed_frequency(seeds_df):
    """Compute the frequency (absolute and %) of each seed compound across genomes."""
    seed_freq = seeds_df.sum(axis=0).sort_values(ascending=False)
    seed_freq_df = seed_freq.reset_index()
    seed_freq_df.columns = ["ModelSEED_ID", "seed_count"]

    total_genomes = seeds_df.shape[0]
    seed_freq_df["frequency_pct"] = (seed_freq_df["seed_count"] / total_genomes) * 100

    return seed_freq_df


def compute_coverage(seeds_df, non_seeds_df):
    """Compute genome coverage (% of genomes) for seeds and non-seeds per compound."""
    total_genomes = seeds_df.shape[0]
    seed_freq = (seeds_df.sum(axis=0) / total_genomes) * 100
    non_seed_freq = (non_seeds_df.sum(axis=0) / total_genomes) * 100

    df_cov = pd.DataFrame(
        {"percent_as_seed": seed_freq, "percent_as_non_seed": non_seed_freq}
    )

    df_cov = df_cov.fillna(0).reset_index().rename(columns={"index": "ModelSEED_ID"})

    return df_cov


# ------------------------
# Visualizations
# ------------------------
def plot_seed_histogram(
    seed_freq_df, bins=30, log_scale=True, color="slateblue", save_path=None
):
    """Histogram: Distribution of seed frequencies across compounds."""
    vals = seed_freq_df["seed_count"].values
    vals = vals[vals > 0]

    if log_scale:
        bins = np.logspace(np.log10(vals.min()), np.log10(vals.max()), bins)

    fig = plt.figure(figsize=(8, 5))
    plt.hist(vals, bins=bins, color=color, edgecolor="black", alpha=0.7)
    if log_scale:
        plt.xscale("log")

    xlabel = "Number of Genomes with Compound as Seed"
    if log_scale:
        xlabel += " (log scale)"
    plt.xlabel(xlabel)

    plt.ylabel("Number of Compounds")
    title = "Distribution of Seed Node Frequencies Across Compounds"
    if log_scale:
        title += " (log scale)"
    plt.title(title)

    if save_path:
        fname = "seed_node_freq_log.png" if log_scale else "seed_node_freq.png"
        plt.savefig(save_path / fname, dpi=300, bbox_inches="tight")

    plt.tight_layout()
    plt.close(fig)
    return None


def plot_top_bar(df, column, top_n=20, color="teal", title=None, save_path=None):
    """Plot a vertical barplot of the top-N compounds by frequency."""
    top = df.nlargest(top_n, column)

    # Use compound names if available
    if "compound_name" in top.columns:
        labels = top["compound_name"]
    else:
        labels = top["ModelSEED_ID"]

    fig = plt.figure(figsize=(12, 5))
    ax = sns.barplot(
        x=labels,
        y=top[column],
        color=color,
    )

    for i, v in enumerate(top[column]):
        ax.text(
            i, v + 0.01 * max(top[column]), str(v), ha="center", va="bottom", fontsize=9
        )
    plt.xlabel("Compound")
    plt.ylabel("Number of Genomes with Compound as Seed")
    plt.title(title if title else f"Top {top_n} Most Essential Metabolites")
    plt.xticks(rotation=45, ha="right")  # Rotate x-axis labels for better readability

    if save_path:
        plt.savefig(
            save_path / f"Top_{top_n}_seed_nodes.png", dpi=300, bbox_inches="tight"
        )

    plt.tight_layout()
    plt.close(fig)
    return None


def plot_coverage_histogram(df_cov, bins=40, logy=False, save_path=None):
    """Plot histogram of genome coverage (%) for Seeds and Non-Seeds compounds."""
    # Prepare bins (0–100%)
    if isinstance(bins, int):
        bins = np.linspace(0, 100, bins + 1)

    fig, ax = plt.subplots(figsize=(7, 5))

    # Histogram seeds
    ax.hist(
        df_cov["percent_as_seed"],
        bins=bins,
        alpha=0.6,
        label="Seeds",
        color="steelblue",
        edgecolor="black",
    )
    # Histogram non-seeds
    ax.hist(
        df_cov["percent_as_non_seed"],
        bins=bins,
        alpha=0.5,
        label="Non-Seeds",
        color="orange",
        edgecolor="black",
    )

    # Add mean lines
    mean_seed = df_cov["percent_as_seed"].mean()
    mean_nonseed = df_cov["percent_as_non_seed"].mean()
    ax.axvline(
        mean_seed, color="red", linestyle="--", label=f"Mean seed = {mean_seed:.1f}%"
    )
    ax.axvline(
        mean_nonseed,
        color="blue",
        linestyle="--",
        label=f"Mean non-seed = {mean_nonseed:.1f}%",
    )

    # Labels
    ax.set_xlabel("Percentage of Genomes")
    ylabel = "Number of Compounds"
    if logy:
        ylabel += " (log scale)"
    ax.set_ylabel(ylabel)

    title = (
        f"Distribution of Genome Coverage (%) for Seeds and Non-Seeds\n"
        f"(Total compounds = {len(df_cov):,})"
    )
    if logy:
        title += " (log scale)"
    ax.set_title(title)

    if logy:
        ax.set_yscale("log")

    ax.legend(loc="upper right", frameon=False, bbox_to_anchor=(1, 1), borderaxespad=0)

    plt.tight_layout()

    # Save
    if save_path:
        fname = (
            "seed_non_seed_genome_coverage_log.png"
            if logy
            else "seed_non_seed_genome_coverage.png"
        )
        out_file = save_path / fname
        plt.savefig(out_file, dpi=300, bbox_inches="tight")
        print(f"Saved plot to {out_file}")

    plt.close(fig)
    return None


# ------------------------
# Main
# ------------------------
def main():
    # Load data
    seeds_df = load_data(SEEDS_PICKLE, filetype="pickle")
    non_seeds_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")

    compounds_map = load_data(COMPOUND_SUMMARY_TSV, filetype="tsv")

    # Compute seed frequencies
    seed_freq_df = compute_seed_frequency(seeds_df)

    # Print quick summary
    print(seed_freq_df[["seed_count", "frequency_pct"]].describe())

    # Ensure output directory exists
    Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

    # Save results
    out_path = Path(OUTPUT_DIR) / "seed_node_frequency.tsv"
    seed_freq_df.to_csv(out_path, sep="\t", index=False)
    print(f"Saved seed node frequency table to {out_path}")

    # Ensure output directory exists
    SEED_NODE_FREQUENCY_PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    save_path = Path(SEED_NODE_FREQUENCY_PLOTS_DIR)

    # Seed histograms (linear + log)
    plot_seed_histogram(seed_freq_df, bins=30, log_scale=False, save_path=save_path)
    plot_seed_histogram(seed_freq_df, bins=30, log_scale=True, save_path=save_path)

    # Merge with compound names for better visualization (optional)
    seed_freq_df = seed_freq_df.merge(
        compounds_map[["ModelSEED_ID", "compound_name"]], on="ModelSEED_ID", how="left"
    )

    # Top-N barplot
    plot_top_bar(seed_freq_df, column="seed_count", top_n=20, save_path=save_path)

    # Coverage histograms (linear + log)
    df_cov = compute_coverage(seeds_df, non_seeds_df)
    df_cov.to_csv(
        Path(OUTPUT_DIR) / "genome_coverage_per_compound.tsv", sep="\t", index=False
    )

    # plot_coverage_histogram(df_cov, bins=40, logy=False, save_path=save_path)
    plot_coverage_histogram(df_cov, bins=40, logy=True, save_path=save_path)


if __name__ == "__main__":
    main()
