"""
Biome-level Metabolic Potential Analysis

This script compares genome-level metabolic potential across soil, marine and freshwater biomes. 
Metabolic potential is quantified by total seeds, total non-seeds, their ratio and values normalized 
per megabase (Mbp). The analysis includes:

1. Statistical tests:
    - Kruskal-Wallis tests to detect biome-level heterogeneity.
    - Dunn's post-hoc tests with FDR correction for pairwise contrasts.
    - Effect sizes (Cliff's delta, median differences) to assess the magnitude of contrasts independently of sample size.

2. Visualizations:
    - Boxplots with significance annotations for each metric across biomes.
    - Stacked barplots by phylum:
        *Mean per-genome metabolic potential (trait-centric).
        *Relative contribution of phyla within each biome (contribution-centric).
    - Heatmap of mean seed/non-seed ratio per phylum per biome.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.stats import kruskal
import scikit_posthocs as sp
import itertools

from config import(
    COMPACT_METADATA_WITH_BIOME_TSV,    # input
    METABOLIC_POTENTIAL_TSV,    # input
    OUTPUT_DIR,
    BIOMES_PLOTS_DIR,
    )
from utils import load_data, split_and_clean_taxonomy, TAXON_PLURALS, LABEL_MAP


# ----------------------------
# Helper Cliff's delta
# ----------------------------
def cliffs_delta(x, y):
    """
    Compute Cliff's delta effect size between two groups.
    0 = no effect
    +/- 1 = complete separation
    """
    nx, ny = len(x), len(y)
    greater = sum(i > j for i in x for j in y)
    less = sum(i < j for i in x for j in y)
    delta = (greater - less) / (nx * ny)
    return delta


# ----------------------------
# Effect sizes and median diffs
# ----------------------------
def compute_effects(df, metric, group_col="main_biome"):
    """Compute pairwise Cliff's delta and median difference for metric across groups."""
    results = []
    groups = df[group_col].dropna().unique()

    for g1, g2 in itertools.combinations(groups, 2):
        x = df.loc[df[group_col] == g1, metric].dropna()
        y = df.loc[df[group_col] == g2, metric].dropna()

        delta = cliffs_delta(x, y)
        med_diff = np.median(x) - np.median(y)

        results.append({
            "Metric": metric,
            "Group1": g1,
            "Group2": g2,
            "Median_Group1": np.median(x),
            "Median_Group2": np.median(y),
            "Median_Diff": med_diff,
            "Cliffs_Delta": delta
        })

    return pd.DataFrame(results)


# -------------------------------
# Boxplot by biome
# -------------------------------
# Annotation
def add_significance_annotations(
    ax, order, posthoc_df, gap=0.06, headroom=0.25, hide_ns=True, fontsize=10
):
    """
    Draw significance brackets for all pairwise comparisons.
    Works with linear or log y-scale.
    """
    ymin, ymax = ax.get_ylim()
    combos = list(itertools.combinations(order, 2))

    # Add headroom so brackets never clip
    if ax.get_yscale() == "log":
        pad = (ymax / ymin) ** headroom
        ax.set_ylim(ymin, ymax * pad)
        y_base = ymax * (1 + gap)
    else:
        pad = (ymax - ymin) * headroom
        ax.set_ylim(ymin, ymax + pad)
        y_base = ymax + (ymax - ymin) * gap

    def stars_from_p(p):
        return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns"

    for i, (g1, g2) in enumerate(combos):
        if g1 not in posthoc_df.index or g2 not in posthoc_df.columns:
            continue
        pval = float(posthoc_df.loc[g1, g2])
        s = stars_from_p(pval)
        if hide_ns and s == "ns":
            continue

        # Vertical level for this bracket (empirical spacing)
        if ax.get_yscale() == "log":
            y = y_base * (1.12 ** i)    
            y_hi = y * 1.06
        else:
            step = (ymax - ymin) * gap
            y = y_base + i * step
            y_hi = y + 0.02 * (ymax - ymin)

        x1, x2 = order.index(g1), order.index(g2)
        ax.plot([x1, x1, x2, x2], [y, y_hi, y_hi, y], lw=1.2, c="black")
        ax.text((x1 + x2) / 2, y_hi, s, ha="center", va="bottom", fontsize=fontsize)


# Boxplot
def plot_boxplot_by_biome(
    df,
    column,
    group_col="main_biome",
    palette="Set2",
    title=None,
    ylabel=None,
    save_dir=None,
    log_scale=False,
    posthoc_df=None,
):
    """Plot boxplot of a metric across biomes (e.g. Soil, Marine, Freshwater)."""
    fig, ax = plt.subplots(figsize=(6, 4))

    if log_scale:
        ax.set_yscale("log")
        ylabel = ylabel or f"Log10({column})"
    else:
        ylabel = ylabel or LABEL_MAP.get(column, column)

    ax.set_title(title if title else f"{column} by Biome")
    ax.set_ylabel(ylabel)

    # x-label from LABEL_MAP
    xlabel = LABEL_MAP.get(group_col, group_col.capitalize())
    ax.set_xlabel(xlabel)
    
    # Enforce consistent order
    order = df[group_col].dropna().unique().tolist()

    sns.boxplot(data=df, x=group_col, y=column, palette=palette, ax=ax, order=order)

    # Sample sizes on x-ticks
    counts = df[group_col].value_counts()
    ax.set_xticks(range(len(order)))
    labels = [f"{b}\n(n={counts.get(b,0)})" for b in order]
    ax.set_xticklabels(labels)

    # Optional significance stars
    if posthoc_df is not None:
        add_significance_annotations(ax, order, posthoc_df, hide_ns=True)

    if save_dir:
        out_path = save_dir / f"boxplot_{column}_by_biome.png"
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"Saved plot: {out_path}")

    plt.tight_layout()
    plt.close(fig)
    return None


# -------------------------------
# Stacked barplot by taxon
# -------------------------------
def plot_stacked_by_taxon(
    df, 
    value_col, 
    taxon="phylum", 
    top_n=10, 
    save_dir=None, 
    log_scale=False,
    mode="mean"  # "mean" or "normalized"
):
    """
    Plot stacked barplot for metabolic potential by taxon across biomes.

    mode:
        - "mean": average value per genome per taxon per biome (trait-centric)
        - "normalized": relative contribution of taxa to each biome (contribution-centric)
    """
    # Extract taxonomy
    taxa = split_and_clean_taxonomy(df, "gtdb_taxonomy")
    df[taxon] = taxa[taxon]

    # Group by biome x taxon
    if mode == "mean":
        grouped = (
            df.groupby(["main_biome", taxon])[value_col]
            .mean()
            .reset_index()
        )
    else:
        grouped = (
            df.groupby(["main_biome", taxon])[value_col]
            .sum()
            .reset_index()
        )

    # Keep only top-N taxa
    top_taxa = df[taxon].value_counts().head(top_n).index
    grouped = grouped[grouped[taxon].isin(top_taxa)]

    if mode == "normalized":
        grouped[value_col] = grouped.groupby("main_biome")[value_col].transform(lambda x: x / x.sum())

    # Pivot for stacked barplot
    pivot = grouped.pivot(index="main_biome", columns=taxon, values=value_col).fillna(0)
    ax = pivot.plot(kind="bar", stacked=True, figsize=(10, 8), colormap="tab20")
    fig = ax.get_figure()

    # Labels
    xlabel = LABEL_MAP.get("main_biome", "main_biome")
    if mode == "mean":
        ylabel = f"Mean {LABEL_MAP.get(value_col, value_col)} per genome"
    else:
        ylabel = "Relative Contribution"

    title_taxon = TAXON_PLURALS.get(taxon, taxon)

    # Plot
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if mode == "mean":
        ax.set_title(f"Metabolic Potential across Biomes (mean) by Top {top_n} {title_taxon}")
    else:
        ax.set_title(f"Relative Contribution of Top {top_n} {title_taxon} to Metabolic Potential across Biomes")

    if log_scale and mode == "mean":
        ax.set_yscale("log")

    ax.legend(
        bbox_to_anchor=(0.5, -0.25),
        loc="upper center",
        ncol=5,
        title=title_taxon
    )

    plt.tight_layout(rect=[0, 0.05, 1, 1])

    if save_dir:
        out_path = save_dir / f"stacked_{mode}_top{top_n}_{taxon}_{value_col}.png"
        plt.savefig(out_path, dpi=300)
        print(f"Saved plot: {out_path}")

    plt.close(fig)
    return None


# -------------------------------
# Heatmap of mean ratios
# -------------------------------
def plot_heatmap_mean_ratio(df, taxon="phylum", top_n=15, save_dir=None):
    """Heatmap of mean seed/non-seed ratio per phylum per biome."""
    # Extract taxonomy
    taxa = split_and_clean_taxonomy(df, "gtdb_taxonomy")
    df[taxon] = taxa[taxon]

    # Filter top taxa by abundance
    top_taxa = df[taxon].value_counts().head(top_n).index
    sub = df[df[taxon].isin(top_taxa)]

    # Compute mean ratios
    grouped = (
        sub.groupby(["main_biome", taxon])["Ratio"]
        .mean()
        .reset_index()
    )

    # Pivot for heatmap
    pivot = grouped.pivot(index=taxon, columns="main_biome", values="Ratio")

    # Plot
    fig, ax = plt.subplots(figsize=(6, 8))
    sns.heatmap(pivot, annot=True, fmt=".2f", cmap="viridis", cbar_kws={'label': 'Mean Ratio'})
    ax.set_xlabel("Biome")
    ax.set_ylabel(taxon.capitalize())
    ax.set_title(f"Mean Seeds/Non-seeds Ratio per {taxon.capitalize()} (Top {top_n})")

    plt.tight_layout()

    if save_dir:
        out_path = save_dir / f"heatmap_mean_ratio_top{top_n}_{taxon}.png"
        plt.savefig(out_path, dpi=300)
        print(f"Saved heatmap: {out_path}")

    plt.close(fig)
    return None


# -------------------------------
# Save stats
# -------------------------------
def save_stats(result_dict, out_path):
    """Save dictionary of stats (metric to results) into TSV file."""
    df = pd.DataFrame(result_dict).T  # metrics as rows
    df.to_csv(out_path, sep="\t", index=True)
    print(f"Saved stats to {out_path}")


# -------------------------------
# Main
# -------------------------------
def main():
    # Load data
    metabolic_df = load_data(METABOLIC_POTENTIAL_TSV, filetype="tsv")

    # Load metadata with biome assignment
    biomes_df = load_data(COMPACT_METADATA_WITH_BIOME_TSV, filetype="tsv")
    biomes_df = biomes_df[biomes_df["main_biome"].isin(["Soil", "Marine", "Freshwater"])]
    biomes_df = biomes_df[["patric_id", "main_biome"]]

    # Merge metabolic data with biome data
    merged = pd.merge(metabolic_df, biomes_df, on="patric_id", how="inner")

    # Metrics
    metrics = [
        "Total_Seeds", "Total_non_Seeds", "Ratio",
        "Seeds_per_Mbp", "Non_Seeds_per_Mbp",
        ]

    # Output directory
    plot_dir = OUTPUT_DIR / BIOMES_PLOTS_DIR
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Collect stats
    kruskal_results = {}
    dunn_results = {}

    for col in metrics:
        print(f"\n=== {col} by biome ===")

        sub = merged.copy()

        if "per_Mbp" in col:
            sub = sub[sub[col] > 0]

        biomes = sub["main_biome"].dropna().unique().tolist()
        groups = [sub.loc[sub["main_biome"] == b, col].dropna() for b in biomes]
        
        # Remove empty groups
        biomes_nonempty = [b for b, g in zip(biomes, groups) if len(g) > 0]
        groups_nonempty = [g for g in groups if len(g) > 0]

        if len(groups_nonempty) < 2:
            print(f"Skipping {col}: <2 non-empty biomes.")
            continue

        # Kruskal-Wallis
        H, p = kruskal(*groups_nonempty)
        kruskal_results[col] = {"H_stat": H, "p_value": p}
        kruskal_results[col].update({b: len(g) for b, g in zip(biomes_nonempty, groups_nonempty)})
        print(f"Kruskal-Wallis H={H:.2f}, p={p:.2e}")

        # Dunn's post-hoc
        posthoc = None
        if len(groups_nonempty) >= 3:
            posthoc = sp.posthoc_dunn(groups_nonempty, p_adjust="fdr_bh")
            posthoc.index = posthoc.columns = biomes_nonempty
            dunn_results[col] = posthoc
            print(posthoc)
        else:
            print("Dunn skipped (only 2 groups).")

        # Plot for each metric
        plot_boxplot_by_biome(
            sub,
            column=col,
            group_col="main_biome",
            palette="Set2",
            save_dir=plot_dir,
            log_scale=("per_Mbp" in col),
            posthoc_df=posthoc,
        )

    # Save stats
    save_stats(kruskal_results, OUTPUT_DIR / "kruskal_results_biomes.tsv")

    # Save all Dunn’s results into one file
    with open(OUTPUT_DIR / "dunn_posthoc_biomes.tsv", "w") as f:
        for col, df in dunn_results.items():
            f.write(f"# Dunn's post-hoc {col}\n")
            df.to_csv(f, sep="\t")
            f.write("\n")
    print("Saved all Dunn’s results to dunn_posthoc_biomes.tsv")

    # Collect effect sizes for all metrics
    effect_dfs = []
    for col in metrics:
        eff = compute_effects(merged, metric=col, group_col="main_biome")
        eff["Metric"] = col
        effect_dfs.append(eff)

    # Save all effect sizes into one file
    with open(OUTPUT_DIR / "effect_sizes_biomes.tsv", "w") as f:
        for eff in effect_dfs:
            metric = eff["Metric"].iloc[0]
            f.write(f"# Effect sizes {metric}\n")
            eff.to_csv(f, sep="\t", index=False)
            f.write("\n")
    print("Saved all effect sizes to effect_sizes_biomes.tsv")

    # Extra plots
    # Stacked barplot (mean)
    plot_stacked_by_taxon(
        merged,
        value_col="Ratio",
        taxon="phylum",
        top_n=15,
        save_dir=plot_dir,
        mode="mean"
    )

    # Stacked barplot (normalized)
    plot_stacked_by_taxon(
        merged,
        value_col="Ratio",
        taxon="phylum",
        top_n=15,
        save_dir=plot_dir,
        mode="normalized"
    )

    # Heatmap of mean ratios
    plot_heatmap_mean_ratio(
        merged,
        taxon="phylum",
        top_n=15,
        save_dir=plot_dir,
    )


if __name__ == "__main__":
    main()
