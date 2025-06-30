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

def report_counts(data):
    """
    Print summary statistics for taxonomic ranks.
    """
    print(f"\nTotal genomes: {len(data)}\n")
    for level in ["domain", "phylum", "class", "order", "family", "genus", "species"]:
        print(f"Unique {level}s: {data[level].nunique()}")


def report_missing_ranks(data):
    """
    Check for empty taxonomic ranks.
    """
    print("\nMissing entries per taxonomic rank:")
    print(data[["domain", "phylum", "class", "order", "family", "genus", "species"]].isna().sum())


def plot_taxonomic_summary(data):
    """
    Generate plots summarizing taxonomy.
    """
    sns.set(style="whitegrid")

    # Bar plot: Unique taxa per rank
    taxa_counts = {
        level.capitalize(): data[level].nunique()
        for level in ["domain", "phylum", "class", "order", "family", "genus", "species"]
    }

    df_taxa = pd.DataFrame({
        "Rank": list(taxa_counts.keys()),
        "Count": list(taxa_counts.values())
    })

    plt.figure(figsize=(10, 5))
    sns.barplot(data=df_taxa, x="Rank", y="Count", hue="Rank", palette="Blues_d", legend=False)
    plt.ylabel("Unique Taxa")
    plt.title("Unique Taxa per Taxonomic Rank")
    plt.tight_layout()
    plt.show()

    # Pie chart: Genomes per domain
    domain_counts = data["domain"].value_counts()
    plt.figure(figsize=(6, 6))
    plt.pie(domain_counts, labels=domain_counts.index,
            autopct="%1.1f%%", startangle=90,
            colors=sns.color_palette("pastel"))
    plt.title("Genomes per Domain")
    plt.axis("equal")
    plt.show()

    # Bar plot: Top 10 phyla
    top_phyla = data["phylum"].value_counts().nlargest(10).reset_index()
    top_phyla.columns = ["Phylum", "Count"]

    plt.figure(figsize=(10, 6))
    sns.barplot(data=top_phyla, x="Count", y="Phylum", hue="Phylum", palette="viridis", legend=False)
    plt.xlabel("Number of Genomes")
    plt.ylabel("Phylum")
    plt.title("Top 10 Most Abundant Phyla")
    plt.tight_layout()
    plt.show()

    # Bar plot: Top 10 genera
    top_genera = data["genus"].value_counts().nlargest(10).reset_index()
    top_genera.columns = ["Genus", "Count"]

    plt.figure(figsize=(10, 6))
    sns.barplot(data=top_genera, x="Count", y="Genus", hue="Genus", palette="magma", legend=False)
    plt.xlabel("Number of Genomes")
    plt.ylabel("Genus")
    plt.title("Top 10 Most Abundant Genera")
    plt.tight_layout()
    plt.show()


def main():
    # Load metadata
    df = load_data(COMPACT_METADATA)
    df = parse_taxonomy(df)
    report_missing_ranks(df)
    report_counts(df)
    plot_taxonomic_summary(df)

if __name__ == "__main__":
    main()
