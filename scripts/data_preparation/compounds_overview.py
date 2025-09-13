"""
Compound Overview and Visualization

This script provides an overview of the metabolic compounds available for downstream analysis.
Starting from two input datasets (Seed and non-Seed binary matrices).

It generates summary visualizations:

- Venn diagram showing compound overlap and exclusivity of compounds between Seed and non-Seed sets
- Histogram of distribution of KEGG modules and pathways per compound
- Bar plot of top KEGG pathways and modules most compounds participate in
- Bar plot of brite ontologies of compounds
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from collections import Counter

from config import (
    SEEDS_PICKLE,  # input
    NON_SEEDS_PICKLE,  # input
    COMPOUND_SUMMARY_TSV,  # input
    OUTPUT_DIR,
    EXPLORATORY_PLOTS_DIR,
)
from utils import load_data, get_compound_sets


# ------------------------
# Explore compounds dataset
# ------------------------
def count_unique(cpd_df, column, return_list=False):
    """Count and display the number of unique KEGG entities (modules, reactions, pathways)
    across all compounds."""
    series = cpd_df[column].dropna().astype(str)

    # Special cleaning for pathways (various delimiters and weird strings)
    if column == "KEGG_pathways":
        series = (
            series.str.replace(";", "|", regex=False)
            .str.replace(",", "|", regex=False)
            .str.replace(" ", "", regex=False)
            .str.split("|")
            .explode()
            .str.strip()
        )
        series = series[~series.isin(["", "nan", "NaN", "None"])]
    else:
        # Modules/Reactions assumed comma-delimited
        series = series.str.split(",").explode().str.strip()
        series = series[series != ""]

    unique_ids = series.unique()
    print(f"Total unique {column}: {len(unique_ids)}\n")

    if return_list:
        print(", ".join(sorted(unique_ids)))

    return unique_ids


# ------------------------
# Summary Statistics
# ------------------------
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
    print(
        f"Distribution (min/25%/50%/75%/max):\n{cpd_df['KEGG_pathway_count'].describe().round(2)}\n"
    )


# ------------------------
# Plot compounds dataset
# ------------------------
def plot_venn(seed_compounds, nonseed_compounds, out_path):
    """Plot Venn diagram showing compound set overlap and unique counts."""
    only_seed = seed_compounds - nonseed_compounds
    only_nonseed = nonseed_compounds - seed_compounds
    common = seed_compounds & nonseed_compounds
    total_unique = len(seed_compounds | nonseed_compounds)

    counts = {
        "Seed Only": len(only_seed),
        "Non-Seed Only": len(only_nonseed),
        "Both": len(common),
    }

    # Plot Venn diagram
    fig = plt.figure(figsize=(7, 5))
    venn2([seed_compounds, nonseed_compounds], set_labels=("Seed", "Non-Seed"))
    plt.title(
        f"Shared and Exclusive Metabolic Compounds between Seed and Non-Seed Set\n"
        f"in the Study Dataset (N = {total_unique:} Compounds)"
    )

    plt.tight_layout()
    plt.savefig(
        out_path / "venn_seed_non_seed_overlap.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)
    return only_seed, only_nonseed, common


# ------------------------
# Histograms
# ------------------------
def plot_kegg_modules_per_compound(cpd_df, out_path):
    """Plot a histogram of the number of KEGG modules associated with each compound."""
    n_compounds = len(cpd_df)
    ccpd_df = cpd_df.copy()
    ccpd_df["KEGG_module_count"] = (
        ccpd_df["KEGG_modules"]
        .fillna("")
        .astype(str)
        .str.split(",")
        .apply(lambda x: len([i for i in x if i.strip()]))
    )

    bins = range(0, ccpd_df["KEGG_module_count"].max() + 2)

    fig = plt.figure(figsize=(8, 4))
    ax = ccpd_df["KEGG_module_count"].hist(
        bins=bins, color="slateblue", edgecolor="black"
    )
    plt.xlabel("Number of KEGG Modules per Compound")
    plt.ylabel("Number of Compounds (log scale)")
    plt.title(
        f"Frequency distribution of metabolic compounds\n across KEGG modules under study (N = {n_compounds} compounds)"
    )
    plt.yscale("log")
    plt.tight_layout()

    # Annotate bars with counts
    for p in ax.patches:
        height = p.get_height()
        if height > 0:
            ax.annotate(
                f"{int(height)}",
                (p.get_x() + p.get_width() / 2, height),
                ha="center",
                va="bottom",
                fontsize=8,
                rotation=0,
            )
    plt.savefig(out_path / "modules_per_cpds.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return None


def plot_kegg_pathways_per_compound(cpd_df, out_path):
    """Plot a histogram of the number of KEGG pathways associated with each compound.
    Shows log y-scale to reveal outliers."""
    n_compounds = len(cpd_df)
    ccpd_df = cpd_df.copy()
    ccpd_df["KEGG_pathway_count"] = (
        ccpd_df["KEGG_pathways"]
        .fillna("")
        .astype(str)
        .str.split(",")
        .apply(lambda x: len([i for i in x if i.strip()]))
    )

    bins = range(0, ccpd_df["KEGG_pathway_count"].max() + 2)

    fig = plt.figure(figsize=(8, 6))
    ax = ccpd_df["KEGG_pathway_count"].hist(bins=bins, color="teal", edgecolor="black")
    plt.xlabel("Number of KEGG Pathways per Compound")
    plt.ylabel("Number of Compounds (log scale)")
    plt.title(
        f"Frequency distribution of metabolic compounds\n across KEGG pathways under study (N = {n_compounds} compounds)"
    )
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
                ha="center",
                va="bottom",
                fontsize=6,
                color="black",
                rotation=90,
            )
        prev_height = height
    plt.savefig(out_path / "pathways_per_cpds.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return None


# ------------------------
# Bar plots
# ------------------------
def get_top_pathways(summary_df, topN=10):
    """Identify top pathways of compounds of study."""
    sub = summary_df.copy()
    path_lists = (
        sub["KEGG_pathways"]
        .dropna()
        .astype(str)
        .str.replace(" ", "", regex=False)
        .str.replace(";", "|", regex=False)
        .str.replace(",", "|", regex=False)
        .str.split("|")
    )
    paths_flat = [p for lst in path_lists for p in lst if p and p.lower() != "nan"]
    counter = Counter(paths_flat)
    return counter.most_common(topN)


def plot_top_pathways(top_pathways, out_path):
    """Barplot to present top pathways of compounds in the study."""
    if not top_pathways:
        print("No pathways found")
        return
    pw_names, pw_counts = zip(*top_pathways)
    fig = plt.figure(figsize=(8, 4))
    ax = plt.barh(pw_names, pw_counts, color="steelblue")

    # Add labels (counts)
    for i, v in enumerate(pw_counts):
        for rect, count in zip(ax.patches, pw_counts):
            plt.text(rect.get_width() + max(pw_counts) * 0.01, rect.get_y() + rect.get_height()/2, str(count), va='center', ha='left')
        plt.text(
            v + max(pw_counts) * 0.01,
            i,
            str(v),
            va="center",
            ha="left",
            fontsize=9
        )

    plt.xlabel("Number of Compounds")
    plt.title("Top KEGG Pathways in the Study Dataset")
    plt.tight_layout()
    plt.savefig(out_path / "top_pathways.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return None


def get_top_modules(summary_df, topN=10):
    """Identify top modules of compounds of study."""
    mod_lists = (
        summary_df["KEGG_modules"]
        .dropna()
        .astype(str)
        .str.replace(" ", "", regex=False)
        .str.split(",")
    )
    mods_flat = [
        m.strip() for lst in mod_lists for m in lst if m and m.strip().lower() != "nan"
    ]
    counter = Counter(mods_flat)
    return counter.most_common(topN)


def plot_top_modules(top_modules, out_path):
    """Barplot to present top modules of compounds in the study."""
    if not top_modules:
        print("No modules found")
        return
    m_names, m_counts = zip(*top_modules)
    fig = plt.figure(figsize=(8, 4))
    ax = plt.barh(m_names, m_counts, color="darkorange")

    # Add labels (counts)
    for i, v in enumerate(m_counts):
        for rect, count in zip(ax.patches, m_counts):
            plt.text(rect.get_width() + max(m_counts) * 0.01, rect.get_y() + rect.get_height()/2, str(count), va='center', ha='left')
        plt.text(
            v + max(m_counts) * 0.01,
            i,
            str(v),
            va="center",
            ha="left",
            fontsize=9
        )

    plt.xlabel("Number of Compounds")
    plt.title("Top KEGG Modules in the Study Dataset")
    plt.tight_layout()
    plt.savefig(out_path / "top_modules.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return None


# ------------------------
# Bar plot
# ------------------------
def plot_brite(cpd_df, out_path, topN=15, show_percent=True):
    """
    Barplot of BRITE ontology (ontology_primary) for compounds.
    topN: number of categories to keep and groups the rest as "Other (rare)".
    """
    if "ontology_primary" not in cpd_df.columns:
        print("Missing column 'ontology_primary' στο compounds_summary.tsv")
        return

    counts = (
        cpd_df["ontology_primary"].fillna("Unclassified").astype(str).value_counts()
    )

    # Keep topN categories
    top = counts.head(topN)
    if len(counts) > topN:
        other_sum = counts.iloc[topN:].sum()
        top.loc["Other (rare)"] = other_sum

    total = top.sum()

    # Plot
    fig = plt.figure(figsize=(10, 6))
    ax = top.sort_values().plot(kind="barh", color="steelblue", edgecolor="black")
    plt.xlabel("Number of compounds")
    plt.ylabel("BRITE ontology")
    plt.title("BRITE ontology distribution across compounds in the Study dataset")

    # Add labels (counts + %)
    for i, v in enumerate(top.sort_values().values):
        if v > 0:
            if show_percent:
                pct = 100 * v / total
                label = f"{v} ({pct:.1f}%)"
            else:
                label = f"{v}"
            ax.text(v + max(top) * 0.01, i, label, va="center", fontsize=9)

    plt.tight_layout()
    plt.savefig(
        Path(out_path) / "brite_ontology_distribution.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)
    return None


# ------------------------
# Main
# ------------------------
def main():
    # Load inputs
    seed_df = load_data(SEEDS_PICKLE, filetype="pickle")
    nonseed_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")
    cpd_df = load_data(COMPOUND_SUMMARY_TSV, filetype="tsv")

    # Summary stats
    compound_summary_stats(cpd_df)

    outfile = Path(OUTPUT_DIR) / "unique_kegg_reactions.tsv"
    outfile.parent.mkdir(parents=True, exist_ok=True)
    out_path = Path(EXPLORATORY_PLOTS_DIR)
    out_path.mkdir(parents=True, exist_ok=True)

    # Unique KEGG entities
    count_unique(cpd_df, "KEGG_modules")
    count_unique(cpd_df, "KEGG_pathways")
    reactions = count_unique(cpd_df, "KEGG_reactions", return_list=False)
    unique_reactions = pd.DataFrame({"KEGG_reaction_ID": sorted(reactions)})
    unique_reactions.to_csv(outfile, sep="\t", index=False)
    print(f"Saved {len(unique_reactions)} unique reactions to {outfile}")

    # Venn
    seed_compounds, nonseed_compounds = get_compound_sets(seed_df, nonseed_df)
    plot_venn(seed_compounds, nonseed_compounds, out_path)

    # Histograms
    plot_kegg_modules_per_compound(cpd_df, out_path)
    plot_kegg_pathways_per_compound(cpd_df, out_path)

    # Top pathways
    top_pathways = get_top_pathways(cpd_df, topN=10)
    plot_top_pathways(top_pathways, out_path)

    # Top modules
    top_modules = get_top_modules(cpd_df, topN=10)
    plot_top_modules(top_modules, out_path)

    # BRITE barplot
    plot_brite(cpd_df, out_path, topN=15, show_percent=True)


if __name__ == "__main__":
    main()
