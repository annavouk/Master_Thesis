"""
Compound Overview and Visualization

This script analyzes metabolic compounds across Seed and Non-Seed datasets
and generates summary visualizations:

- Bar chart and Venn diagram showing compound overlap between Seed and Non-Seed sets
- Distribution of KEGG modules and pathways per compound
- Heatmap of compound participation across KEGG pathways

Data is read from preprocessed summary and matrix files.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from utils import load_data
from config import (
    SEEDS_PICKLE,
    NON_SEEDS_PICKLE,
    COMPOUND_SUMMARY_TSV,
    )


def get_compound_sets(seed_df, nonseed_df):
    """
    Extract sets of compounds from each matrix.

    Parameters:
        seed_df (pd.DataFrame): Seed binary matrix.
        nonseed_df (pd.DataFrame): Non-seed binary matrix.

    Returns:
        tuple: (set of seed compounds, set of non-seed compounds)
    """
    return set(seed_df.columns), set(nonseed_df.columns)


def plot_bar_and_venn(seed_compounds, nonseed_compounds):
    """
    Plot bar chart and Venn diagram showing compound set overlap.

    Parameters:
        seed_compounds (set): Set of compounds in SeedSet.
        nonseed_compounds (set): Set of compounds in NonSeedSet.
    """
    only_seed = seed_compounds - nonseed_compounds
    only_nonseed = nonseed_compounds - seed_compounds
    common = seed_compounds & nonseed_compounds

    counts = {
        "Only in SeedSet": len(only_seed),
        "Only in NonSeedSet": len(only_nonseed),
        "In Both": len(common),
    }

    plt.figure(figsize=(6, 4))
    plt.bar(counts.keys(), counts.values(), color=["skyblue", "lightcoral", "mediumseagreen"])
    plt.ylabel("Number of Compounds")
    plt.title("Compound Distribution Across Sets")
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(5, 5))
    venn2([seed_compounds, nonseed_compounds], set_labels=("SeedSet", "NonSeedSet"))
    plt.title("Compound Overlap")
    plt.tight_layout()
    plt.show()


def count_unique_modules(cpd_df):
    """
    Count and display the number of unique KEGG modules across all compounds.

    Parameters:
        cpd_df (pd.DataFrame): Compound summary dataframe.
    """
    # Drop missing values and split on comma
    all_modules = cpd_df['KEGG_modules'].dropna().str.split(",").explode().str.strip()

    # Drop empty strings and count uniques
    unique_modules = all_modules[all_modules != ""].unique()
    print(f"Total unique KEGG modules: {len(unique_modules)}\n")
    print("Module IDs:")
    print(", ".join(sorted(unique_modules)))


def count_unique_pathways(cpd_df):
    """
    Count and display the number of unique KEGG pathways across all compounds.

    Parameters:
        cpd_df (pd.DataFrame): Compound summary dataframe.
    """
    # Drop missing values and split on comma
    all_pathways = cpd_df['KEGG_pathways'].dropna().str.split(",").explode().str.strip()

    # Drop empty strings and count uniques
    unique_pathways = all_pathways[all_pathways != ""].unique()
    print(f"Total unique KEGG pathways: {len(unique_pathways)}\n")
    print("Pathway IDs:")
    print(", ".join(sorted(unique_pathways)))


def compound_summary_stats(cpd_df):
    cpd_df["KEGG_pathway_count"] = (
        cpd_df["KEGG_pathways"]
        .fillna("")
        .astype(str)
        .str.split(",")
        .apply(lambda x: len([i for i in x if i.strip()]))
    )
    avg_pathways = cpd_df["KEGG_pathway_count"].mean()
    median_pathways = cpd_df["KEGG_pathway_count"].median()
    orphan_count = (cpd_df["KEGG_pathway_count"] == 0).sum()
    orphan_percent = orphan_count / len(cpd_df) * 100

    print("\nSummary Statistics")
    print(f"Total compounds: {len(cpd_df)}")
    print(f"Mean pathways per compound: {avg_pathways:.2f}")
    print(f"Median pathways per compound: {median_pathways:.2f}")
    print(f"Orphan compounds (no pathway): {orphan_count} ({orphan_percent:.1f}%)")
    print(f"Distribution (min/25%/50%/75%/max):\n{cpd_df['KEGG_pathway_count'].describe().round(2)}\n")


def plot_kegg_modules_per_compound(cpd_df):
    """
    Plot the number of KEGG modules associated with each compound.
    """
    cpd_df["KEGG_module_count"] = cpd_df["KEGG_modules"].str.split(",").apply(lambda x: len([i for i in x if i.strip()]))
    plt.figure(figsize=(6, 4))
    cpd_df["KEGG_module_count"].hist(bins=20, color="slateblue", edgecolor="black")
    plt.xlabel("Number of KEGG Modules")
    plt.ylabel("Number of Compounds")
    plt.title("KEGG Modules per Compound")
    plt.tight_layout()
    plt.show()


def plot_kegg_pathways_per_compound(cpd_df):
    """
    Plot the number of KEGG pathways associated with each compound.
    """
    cpd_df["KEGG_pathway_count"] = (
        cpd_df["KEGG_pathways"]
        .fillna("")  # Replace NaNs with empty strings
        .astype(str)  # Ensure all entries are strings
        .str.split(",")
        .apply(lambda x: len([i for i in x if i.strip()]))  # Count non-empty items
    )

    plt.figure(figsize=(6, 4))
    cpd_df["KEGG_pathway_count"].hist(bins=20, color="teal", edgecolor="black")
    plt.xlabel("Number of KEGG Pathways")
    plt.ylabel("Number of Compounds")
    plt.title("KEGG Pathways per Compound")
    plt.tight_layout()
    plt.show()


def main():
    seed_df = load_data(SEEDS_PICKLE, filetype="pickle")
    nonseed_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")

    cpd_df = load_data(COMPOUND_SUMMARY_TSV, filetype="tsv")

    seed_compounds, nonseed_compounds = get_compound_sets(seed_df, nonseed_df)
    plot_bar_and_venn(seed_compounds, nonseed_compounds)

    count_unique_modules(cpd_df)
    count_unique_pathways(cpd_df)
    compound_summary_stats(cpd_df)
    plot_kegg_modules_per_compound(cpd_df)
    plot_kegg_pathways_per_compound(cpd_df)


if __name__ == "__main__":
    main()
