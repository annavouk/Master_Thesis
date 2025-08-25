"""
Metabolic Potential - Pathway Coverage Visualization (Approach 2)

This script visualizes the metabolic functional potential of genomes at the pathway level.
It uses genome-by-pathway coverage results (fraction of pathway compounds present as non-seeds per genome)
to generate summary plots, including:

- Heatmap of pathway coverage for top genomes and pathways
- Scatter plot: Number of fully covered pathways per genome vs. genome size
- Boxplot and violin plot: Distribution of pathway coverage by phylum
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from config import GENOME_PATHWAY_COVERAGE, COMPACT_METADATA
from utils import load_data, split_and_clean_taxonomy 

def preprocess_coverage(coverage_df, metadata, coverage_threshold=0.9):
    """Calculate the number of fully covered pathways per genome and merge with metadata."""
    high_cov = coverage_df[coverage_df["coverage"] >= coverage_threshold]
    n_pathways = high_cov.groupby("patric_id")["KEGG_pathway"].count().reset_index(name="n_full_pathways")
    merged = n_pathways.merge(metadata[["patric_id", "genome_length", "gtdb_taxonomy"]], on="patric_id", how="left")
    # Taxonomy split: Εδώ βάλε ΜΟΝΟ ΕΝΑ helper, όχι δύο.
    taxonomy_df = split_and_clean_taxonomy(merged, "gtdb_taxonomy")
    merged = pd.concat([merged, taxonomy_df], axis=1)
    return merged

def plot_heatmap_coverage(coverage_df, n_pathways=20, n_genomes=40):
    """Plot a heatmap of pathway coverage for the top N genomes and pathways."""

    top_pathways = coverage_df["KEGG_pathway"].value_counts().nlargest(n_pathways).index
    top_genomes = coverage_df["patric_id"].value_counts().nlargest(n_genomes).index
    df_sub = coverage_df[
        (coverage_df["KEGG_pathway"].isin(top_pathways)) &
        (coverage_df["patric_id"].isin(top_genomes))
    ].drop_duplicates(subset=["patric_id", "KEGG_pathway"], keep="first")
    heatmap_data = df_sub.pivot(index="patric_id", columns="KEGG_pathway", values="coverage")
    plt.figure(figsize=(12, 8))
    sns.heatmap(heatmap_data, cmap="YlGnBu", cbar_kws={'label': 'Coverage'})
    plt.title("Pathway Coverage per Genome (Top Genomes & Pathways)")
    plt.xlabel("KEGG Pathway")
    plt.ylabel("Genome (patric_id)")
    plt.tight_layout()
    plt.show()

def plot_scatter_pathways_vs_genomes(merged, phylum_col="phylum"):
    """Plot a scatter plot: Number of fully covered pathways per genome vs genome size."""
    plt.figure(figsize=(8,5))
    sns.scatterplot(data=merged, x="genome_length", y="n_full_pathways", hue=phylum_col, alpha=0.6)
    plt.xlabel("Genome Length (bp)")
    plt.ylabel("# of Fully Covered Pathways")
    plt.title("Full Pathway Coverage vs Genome Size")
    plt.tight_layout()
    plt.show()

def plot_boxplot_pathways_by_phylum(merged, phylum_col="phylum", top_n=7):
    """Boxplot: Distribution of fully covered pathways by phylum."""
    plt.figure(figsize=(10,5))
    top_phyla = merged[phylum_col].value_counts().nlargest(top_n).index
    sns.boxplot(data=merged[merged[phylum_col].isin(top_phyla)], x=phylum_col, y="n_full_pathways")
    plt.xlabel("Phylum")
    plt.ylabel("# of Fully Covered Pathways")
    plt.title("Distribution of Full Pathway Coverage by Phylum")
    plt.tight_layout()
    plt.show()

def plot_violin_pathways_by_phylum(merged, phylum_col="phylum", top_n=7):
    """Violin plot: Distribution of fully covered pathways by phylum."""
    plt.figure(figsize=(12,6))
    top_phyla = merged[phylum_col].value_counts().nlargest(top_n).index
    sns.violinplot(data=merged[merged[phylum_col].isin(top_phyla)], x=phylum_col, y="n_full_pathways", inner="box")
    plt.xlabel("Phylum")
    plt.ylabel("# of Fully Covered Pathways")
    plt.title("Violin Plot: Full Pathway Coverage by Phylum")
    plt.tight_layout()
    plt.show()

def main():
    # Load data
    coverage_df = load_data(GENOME_PATHWAY_COVERAGE, filetype="tsv")
    metadata = load_data(COMPACT_METADATA, filetype="csv")
    
    # Preprocess and merge  
    merged = preprocess_coverage(coverage_df, metadata)

    # Vusualize
    plot_heatmap_coverage(coverage_df)
    plot_scatter_pathways_vs_genomes(merged)
    plot_boxplot_pathways_by_phylum(merged)
    plot_violin_pathways_by_phylum(merged)

if __name__ == "__main__":
    main()
