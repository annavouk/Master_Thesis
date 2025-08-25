"""
Overview of Genome Taxonomy Distribution

This script summarizes the taxonomic distribution of a genome dataset.
It extracts taxonomic ranks from GTDB-formatted taxonomy strings and generates
summary statistics and visualizations, including barplots and pie charts for key ranks.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from utils import load_data, parse_taxonomy
from config import COMPACT_METADATA


# ------------------------
# Report counts
# ------------------------
def report_counts(data):
    """Print the length of each taxonomic rank."""
    print(f"\nTotal genomes: {len(data)}\n")
    for level in ["domain", "phylum", "class", "order", "family", "genus", "species"]:
        print(f"Unique {level}: {data[level].nunique():,}")


def report_missing_ranks(data):
    """Check for empty taxonomic ranks."""
    print("\nMissing entries per taxonomic rank:")
    print(
        data[["domain", "phylum", "class", "order", "family", "genus", "species"]]
        .isna()
        .sum()
    )


# ------------------------
# Plot
# ------------------------
def plot_taxonomic_summary(data):
    """Generate plots summarizing taxonomy."""
    total_genomes = len(data)
    sns.set(style="whitegrid")

    # Bar plot: Unique taxa per rank
    taxa_counts = {
        level.capitalize(): data[level].nunique()
        for level in [
            "domain",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ]
    }

    df_taxa = pd.DataFrame(
        {"Rank": list(taxa_counts.keys()), "Count": list(taxa_counts.values())}
    )

    # Create two bar plots in a single figure
    fig, (ax1, ax2) = plt.subplots(
        2, 1, sharex=True, gridspec_kw={"height_ratios": [1, 3]}, figsize=(10, 7)
    )

    # Plot the same barplot on two axes (for broken y-axis effect)
    sns.barplot(
        data=df_taxa,
        x="Rank",
        y="Count",
        hue="Rank",
        palette="Blues_d",
        legend=False,
        ax=ax1,
    )
    sns.barplot(
        data=df_taxa,
        x="Rank",
        y="Count",
        hue="Rank",
        palette="Blues_d",
        legend=False,
        ax=ax2,
    )

    # Set the limits on y axis
    species_count = df_taxa["Count"].iloc[-1]
    other_max = df_taxa["Count"].iloc[:-1].max()

    ax1.set_ylim(species_count - 3000, species_count + 2000)
    ax2.set_ylim(0, other_max + 1000)

    # Break mark
    ax1.spines["bottom"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax1.tick_params(labeltop=False)
    ax2.xaxis.tick_bottom()

    d = 0.008
    kwargs = dict(transform=ax1.transAxes, color="k", clip_on=False)
    ax1.plot((-d, +d), (-d, +d), **kwargs)
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
    kwargs.update(transform=ax2.transAxes)
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    # Use the y-limit break as the threshold
    break_threshold = species_count - 3000

    # Place value labels above each bar (upper axis)
    bars1 = ax1.patches
    for bar, row in zip(bars1, df_taxa.itertuples()):
        value = row.Count
        if not value or pd.isna(value):
            continue  # Skip bars with zero or missing count
        label = f"{value:,}"  # Add thousands separator
        # Capped offset for label
        offset = min(80, max(0.03 * value, 10))
        x = bar.get_x() + bar.get_width() / 2
        y = bar.get_height()
        if value > break_threshold:
            ax1.text(
                x,
                y + offset,
                label,
                ha="center",
                va="bottom",
                fontsize=12,
                fontweight="bold",
            )

    # Place value labels above each bar (lower axis)
    bars2 = ax2.patches
    for bar, row in zip(bars2, df_taxa.itertuples()):
        value = row.Count
        if not value or pd.isna(value):
            continue  # Skip bars with zero or missing count
        label = f"{value:,}"
        offset = min(80, max(0.03 * value, 10))
        x = bar.get_x() + bar.get_width() / 2
        y = bar.get_height()
        if value <= break_threshold:
            ax2.text(
                x,
                y + offset,
                label,
                ha="center",
                va="bottom",
                fontsize=12,
                fontweight="bold",
            )

    fig.suptitle(
        f"Number of Unique Taxonomic Groups per Rank in the Study Dataset (N = {total_genomes:,} Genomes)",
        y=0.98,
        fontsize=15,
    )
    ax2.set_ylabel("Number of Unique Taxonomic Groups")
    ax2.set_xlabel("Taxonomic Rank")
    ax1.set_ylabel("")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

    # Pie chart: Genomes per domain
    domain_counts = data["domain"].value_counts()
    total = domain_counts.sum()

    plt.figure(figsize=(6, 6))
    plt.pie(
        domain_counts,
        labels=domain_counts.index,
        autopct=lambda pct: f"{int(pct*total/100):,} ({pct:.1f}%)",
        startangle=90,
        colors=sns.color_palette("pastel"),
        textprops={"fontsize": 12},
    )
    plt.title(
        f"Distribution of Genomes Across Domains in the Study Dataset\n(N = {total_genomes:,} Genomes)",
        y=1.05
    )
    plt.axis("equal")
    plt.tight_layout()
    plt.show()

    # Bar plot: Top 10 phyla
    top_phyla = data["phylum"].value_counts().nlargest(10).reset_index()
    top_phyla.columns = ["Phylum", "Count"]

    top_phyla_total = top_phyla["Count"].sum()

    plt.figure(figsize=(10, 6))
    ax = sns.barplot(
        data=top_phyla,
        x="Count",
        y="Phylum",
        hue="Phylum",
        palette="viridis",
        legend=False,
    )
    plt.xlabel("Number of Genomes")
    plt.ylabel("Phylum")
    plt.title(
        f"Top 10 Most Abundant Phyla in the Study Dataset\n"
        f"({top_phyla_total:,} of {total_genomes:,} Genomes)"
    )

    for i, row in top_phyla.iterrows():
        ax.text(
            row["Count"] + max(top_phyla["Count"]) * 0.01,
            i,
            f'{row["Count"]:,}',
            va="center",
            fontsize=12,
        )

    plt.tight_layout()
    plt.show()

    # Bar plot: Top 10 genera
    top_genera = data["genus"].value_counts().nlargest(10).reset_index()
    top_genera.columns = ["Genus", "Count"]

    top_genera_total = top_genera["Count"].sum()
    total_genomes = len(data)

    plt.figure(figsize=(10, 6))
    ax = sns.barplot(
        data=top_genera,
        x="Count",
        y="Genus",
        hue="Genus",
        palette="magma",
        legend=False,
    )
    plt.xlabel("Number of Genomes")
    plt.ylabel("Genus")
    plt.title(
        f"Top 10 Most Abundant Genera in the Study Dataset\n"
        f"({top_genera_total:,} of {total_genomes:,} Genomes)"
    )

    for i, row in top_genera.iterrows():
        ax.text(
            row["Count"] + max(top_genera["Count"]) * 0.01,
            i,
            f'{row["Count"]:,}',
            va="center",
            fontsize=12,
        )

    plt.tight_layout()
    plt.show()


# ------------------------
# Main
# ------------------------
def main():
    # Load metadata
    df = load_data(COMPACT_METADATA)
    df = parse_taxonomy(df)
    report_missing_ranks(df)
    report_counts(df)
    plot_taxonomic_summary(df)


if __name__ == "__main__":
    main()
