"""
Compound Overview and Visualization

This script provides an overview of the metabolic compounds available for downstream analysis.
Starting from two input datasets (Seed and Non-Seed binary matrices), we first assess the overlap and exclusivity of compounds between them.

compound summary tsv file

generates summary visualizations:

- Bar chart and Venn diagram showing compound overlap and exclusivity of compounds between Seed and Non-Seed sets
- Histogram of distribution of KEGG modules and pathways per compound
- Bar plot of top KEGG pathways in which participate the compounds from Seed, non-Seed subset and their overlap.

Data is read from preprocessed summary and matrix files.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from collections import Counter
from utils import load_data
from config import (
    SEEDS_PICKLE,
    NON_SEEDS_PICKLE,
    COMPOUND_SUMMARY_TSV,
    )


def get_compound_sets(seed_df, nonseed_df):
    """
    Extract sets of compounds from each binary matrix.

    Parameters:
        seed_df (pd.DataFrame): Seed binary matrix.
        nonseed_df (pd.DataFrame): Non-seed binary matrix.

    Returns:
        tuple: (set of seed compounds, set of non-seed compounds)
    """
    return set(seed_df.columns), set(nonseed_df.columns)


def plot_bar_and_venn(seed_compounds, nonseed_compounds):
    """
    Plot a bar chart and a Venn diagram showing compound set overlap.

    Parameters:
        seed_compounds (set): Set of compounds in SeedSet.
        nonseed_compounds (set): Set of compounds in NonSeedSet.
    """
    only_seed = seed_compounds - nonseed_compounds
    only_nonseed = nonseed_compounds - seed_compounds
    common = seed_compounds & nonseed_compounds
    total_unique = len(seed_compounds | nonseed_compounds)

    counts = {
        "Seed Only": len(only_seed),
        "Non-Seed Only": len(only_nonseed),
        "Both": len(common),
    }

    # Plot bar chart
    plt.figure(figsize=(8, 6))
    bars = plt.bar(counts.keys(), counts.values(), color=["skyblue", "lightcoral", "mediumseagreen"])
    plt.ylabel("Number of Unique Compounds")
    plt.title(f"Number of Unique Metabolic Compounds Classified as Seed Only, Non-Seed Only, or Both\n"
    f"in the Study Dataset (N = {total_unique:} Compounds)")

    for bar in bars:
        height = bar.get_height()
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            height + max(counts.values())*0.01,   
            f"{int(height):}",
            ha='center', va='bottom', fontsize=12, fontweight='bold'
        )

    plt.tight_layout()
    plt.show()

    # Plot Venn diagram
    plt.figure(figsize=(7, 5))
    venn2([seed_compounds, nonseed_compounds], set_labels=("Seed", "Non-Seed"))
    plt.title(f"Shared and Exclusive Metabolic Compounds between Seed and Non-Seed Set\n"
    f"in the Study Dataset (N = {total_unique:} Compounds)")

    plt.tight_layout()
    plt.show()

    return only_seed, only_nonseed, common

def count_unique_modules(cpd_df):
    """
    Count and display the number of unique KEGG modules across all compounds.

    Parameters:
        cpd_df (pd.DataFrame): Compound summary dataframe with 'KEGG_modules' column.
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
        cpd_df (pd.DataFrame): Compound summary dataframe 'KEGG_pathways' column.
    """
    # Drop missing values and split on comma
    all_pathways = cpd_df['KEGG_pathways'].dropna().str.split(",").explode().str.strip()

    # Drop empty strings and count uniques
    unique_pathways = all_pathways[all_pathways != ""].unique()
    print(f"Total unique KEGG pathways: {len(unique_pathways)}\n")
    print("Pathway IDs:")
    print(", ".join(sorted(unique_pathways)))


def compound_summary_stats(cpd_df):
    """
    Calculates and prints summary statistics for the number of KEGG pathways associated with each compound.

    For each compound, counts how many KEGG pathways it is annotated with,
    then reports:
      - Total number of compounds
      - Mean and median pathways per compound
      - Number (and percentage) of 'orphan' compounds (not annotated in any pathway)
      - Distribution summary (min, 25th percentile, median, 75th percentile, max)
    """
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
    Plot a histogram of the number of KEGG modules associated with each compound.
    """
    n_compounds = len(cpd_df)
    #cpd_df = cpd_df.copy()
    cpd_df["KEGG_module_count"] = cpd_df["KEGG_modules"].fillna("").astype(str).str.split(",").apply(lambda x: len([i for i in x if i.strip()]))
    
    bins = range(1, cpd_df["KEGG_module_count"].max() + 2)
    
    plt.figure(figsize=(8, 4))
    ax = cpd_df["KEGG_module_count"].hist(bins=bins, color="slateblue", edgecolor="black")
    plt.xlabel("Number of KEGG Modules per Compound")
    plt.ylabel("Number of Compounds (log scale)")
    plt.title(f"Distribution of KEGG Module Count per Metabolic Compound\nin the Study Dataset (N = {n_compounds} compounds)")
    plt.yscale('log')
    plt.tight_layout()
   
   # Annotate bars with counts
    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.annotate(f'{int(height)}', (p.get_x() + p.get_width() / 2, height),
                        ha='center', va='bottom', fontsize=8, rotation=0)
    plt.show()

def plot_kegg_pathways_per_compound(cpd_df):
    """
    Plot a histogram of the number of KEGG pathways associated with each compound.
    Shows log y-scale to reveal outliers.
    """
    n_compounds = len(cpd_df)
    #cpd_df = cpd_df.copy()
    cpd_df["KEGG_pathway_count"] = cpd_df["KEGG_pathways"].fillna("").astype(str).str.split(",").apply(lambda x: len([i for i in x if i.strip()]))

    bins = range(0, cpd_df["KEGG_pathway_count"].max() + 2)

    plt.figure(figsize=(8, 6))
    ax = cpd_df["KEGG_pathway_count"].hist(bins=bins, color="teal", edgecolor="black")
    plt.xlabel("Number of KEGG Pathways per Compound")
    plt.ylabel("Number of Compounds (log scale)")
    plt.title(f"KEGG Pathways per Compound\nin the Study Dataset (N = {n_compounds} compounds)")
    plt.yscale("log")
    plt.tight_layout()

    # Annotate bars with counts and add label only to the first occurrence of each height 
    prev_height = None
    for p in ax.patches:
        height = p.get_height()
        offset = 0.05 * height
        if height > 10 and height != prev_height:
            ax.annotate(
                f"{int(height)}",
                (p.get_x() + p.get_width() / 2, height + offset),
                ha="left", va="bottom", fontsize=7, color="black", rotation=90
            )
        prev_height = height
    plt.show()

def get_top_pathways(compounds_set, summary_df, topN=10):
    """
    Calculate the most frequent KEGG pathways among a given set of compounds.

    Parameters:
        compounds_set (set): Set of compound IDs to analyze.
        summary_df (pd.DataFrame): DataFrame containing compound metadata. 
                                   Must include columns 'SEED_ID' and 'KEGG_pathways'.
        topN (int): Number of top pathways to return.

    Returns:
        list: List of pairs (pathway, count) for the topN most frequent pathways.
    """
    sub = summary_df[summary_df['SEED_ID'].isin(compounds_set)]
    # Assume pathways is string, split on ; or | or ,
    path_lists = sub['KEGG_pathways'].dropna().astype(str).str.replace(' ', '').str.replace(';', '|').str.replace(',', '|').str.split('|')
    paths_flat = [p for sublist in path_lists for p in sublist if p and p != 'nan']
    counter = Counter(paths_flat)
    return counter.most_common(topN)

def plot_top_pathways(top_pathways, subset_name, cpd_df, compounds_set):
    """
    Plot a horizontal bar plot of the most frequent KEGG pathways in a given compound subset.

    Parameters:
        top_pathways : list of (pathway, count) pairs.
        subset_name : str
            Name/label for the compound subset (for plot title).
    """
    
    if not top_pathways:
        print(f"No pathways found for {subset_name}")
        return

    # All unique pathways in the dataset
    all_pathways = cpd_df['KEGG_pathways'].dropna().str.split(",").explode().str.strip()
    all_pathways = all_pathways[all_pathways != ""]
    total_unique = len(all_pathways.unique())

    # Unique pathways in the subset
    subset_df = cpd_df[cpd_df['SEED_ID'].isin(compounds_set)]
    subset_pathways = subset_df['KEGG_pathways'].dropna().str.split(",").explode().str.strip()
    subset_pathways = subset_pathways[subset_pathways != ""]
    subset_unique = len(subset_pathways.unique())

    pw_names, pw_counts = zip(*top_pathways)
    plt.figure(figsize=(8, 4))
    plt.barh(pw_names, pw_counts, color='steelblue')
    plt.xlabel('Number of Compounds')
    plt.title(
        f"Top KEGG Pathways by Number of Compounds Functioning as {subset_name}\n"
        f"(Total Pathways in Dataset: {total_unique}, "
        f"Total in {subset_name} Subset: {subset_unique})"
    )
    plt.tight_layout()
    plt.show()


def main():
    seed_df = load_data(SEEDS_PICKLE, filetype="pickle")
    nonseed_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")

    cpd_df = load_data(COMPOUND_SUMMARY_TSV, filetype="tsv")

    seed_compounds, nonseed_compounds = get_compound_sets(seed_df, nonseed_df)
    only_seed, only_nonseed, common = plot_bar_and_venn(seed_compounds, nonseed_compounds)
    #print(f"Compounds identified as seeds only:", only_seed)
    #print(f"Compounds identified as non-seeds only:", only_nonseed)
    #print(f"Compounds identified as both seeds and non-seeds:", common)

    #seed_only_df = cpd_df[cpd_df["SEED_ID"].isin(only_seed)]
    #count_unique_pathways(seed_only_df)

    #count_unique_modules(cpd_df)
    #count_unique_pathways(cpd_df)
    #compound_summary_stats(cpd_df)
    
    plot_kegg_modules_per_compound(cpd_df)
    #module_count_freq = cpd_df["KEGG_module_count"].value_counts().sort_index()
    #print(module_count_freq)
    #cpds_15_30 = cpd_df[cpd_df["KEGG_module_count"].between(15, 30)]
    #print(cpds_15_30[["SEED_ID", "KEGG_modules", "KEGG_module_count", "dataset"]])

    #plot_kegg_pathways_per_compound(cpd_df)
    #pathway_count_freq = cpd_df["KEGG_pathway_count"].value_counts().sort_index()
    #print(pathway_count_freq)


    #for subset_name, subset_set in zip(['Seed Only', 'Non-Seed Only', 'Seed and Non-Seed'], 
    #                                   [only_seed, only_nonseed, common]):
    #    top_pathways = get_top_pathways(subset_set, cpd_df, topN=10)
    #    print(f"\nTop KEGG pathways for {subset_name}:")
    #    for pw, count in top_pathways:
    #        print(f"{pw}: {count}")
    #    plot_top_pathways(top_pathways, subset_name, cpd_df, subset_set)


if __name__ == "__main__":
    main()
