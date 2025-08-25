"""
Analyze the frequency of seed nodes (essential metabolites) across genomes.

Generate various plots (histogram,) to visualize the distribution of seed node frequencies.
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from config import SEEDS_PICKLE
from utils import load_data, plot_histogram


# ------------------------
# Compute seed node frequency 
# ------------------------
def compute_seed_frequency(seeds_df):
    """Compute the frequency of each seed compound across genomes."""
    seed_freq = seeds_df.sum(axis=0).sort_values(ascending=False)
    seed_freq_df = seed_freq.reset_index()
    seed_freq_df.columns = ['compound', 'seed_count']
    return seed_freq_df


# ------------------------
# Visualizations
# ------------------------
def plot_seed_histogram(seed_freq_df, bins=30, log_scale=True, color="slateblue"):
    """Histogram: Distribution of seed frequencies across compounds."""
    vals = seed_freq_df['seed_count'].values
    vals = vals[vals > 0]

    if log_scale:
        bins = np.logspace(np.log10(vals.min()), np.log10(vals.max()), bins)

    plt.figure(figsize=(8,5))
    plt.hist(vals, bins=bins, color=color, edgecolor="black", alpha=0.7)
    if log_scale:
        plt.xscale("log")

    plt.xlabel("Number of Genomes with Compound as Seed")
    plt.ylabel("Number of Compounds")
    plt.title("Distribution of Seed Node Frequencies Across Compounds")
    plt.tight_layout()
    plt.show()


def plot_top_bar(df, column, top_n=20, color="teal", title=None):
    """Plot a vertical barplot of the top-N compounds by frequency."""
    top = df.nlargest(top_n, column)
    plt.figure(figsize=(12, 5))
    ax = sns.barplot(
        x=top['compound'],      # compound names στον x-άξονα
        y=top[column],          # counts στον y-άξονα
        color=color,
    )
    for i, v in enumerate(top[column]):
        ax.text(i, v + 0.01 * max(top[column]), str(v),
                ha='center', va='bottom', fontsize=9)
    plt.xlabel("Compound (ModelSEED ID)")
    plt.ylabel("Number of Genomes with Compound as Seed")
    plt.title(title if title else f"Top {top_n} Most Essential Metabolites")
    plt.xticks(rotation=45, ha="right")  # για να μη στριμώχνονται τα labels
    plt.tight_layout()
    plt.show()


# ------------------------
# Main
# ------------------------
def main():
    # Load data
    seeds_df = load_data(SEEDS_PICKLE, filetype='pickle')

    # Compute seed frequencies
    seed_freq_df = compute_seed_frequency(seeds_df)
    print(seed_freq_df['seed_count'].describe())

    # Histogram of seed frequencies
    plot_seed_histogram(seed_freq_df, bins=30, log_scale=True)

    # Barplot of top-N seed compounds
    plot_top_bar(seed_freq_df, column="seed_count", top_n=20)


if __name__ == "__main__":
    main()
