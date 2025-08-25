"""
Compare genome-level metabolic potential across biomes (Soil, Marine, Freshwater).

Generates boxplots and barplots.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal
import scikit_posthocs as sp
import itertools

from config import COMPACT_METADATA_WITH_BIOME, METABOLIC_POTENTIAL_1, OUTPUT_DIR
from utils import load_data, split_and_clean_taxonomy, TAXON_PLURALS, LABEL_MAP


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

    # add headroom so brackets never clip
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

        # vertical level for this bracket
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
    """Plot boxplot of a metric across biomes (e.g., Soil, Marine, Freshwater)."""
    fig, ax = plt.subplots(figsize=(6, 4))

    if log_scale:
        ax.set_yscale("log")
        ylabel = ylabel or f"Log10({column})"
    else:
        ylabel = ylabel or LABEL_MAP.get(column, column)

    ax.set_title(title if title else f"{column} by Biome")
    ax.set_ylabel(ylabel)

    # X-label from LABEL_MAP
    xlabel = LABEL_MAP.get(group_col, group_col.capitalize())
    ax.set_xlabel(xlabel)

    # Annotate sample sizes
    counts = df[group_col].value_counts()
    labels = [f"{b}\n(n={counts[b]})" for b in df[group_col].unique()]
    ax.set_xticklabels(labels)
    
    # Enforce consistent order
    order = ["Soil", "Marine", "Freshwater"]

    sns.boxplot(data=df, x=group_col, y=column, palette=palette, ax=ax, order=order)

    # Sample sizes on x-ticks
    counts = df[group_col].value_counts()
    labels = [f"{b}\n(n={counts.get(b,0)})" for b in order]
    ax.set_xticklabels(labels)

    # Optional significance stars
    if posthoc_df is not None:
        add_significance_annotations(ax, order, posthoc_df, hide_ns=True)

    plt.tight_layout()

    if save_dir:
        out_path = save_dir / f"boxplot_{column}_by_biome.png"
        plt.savefig(out_path, dpi=300)
        print(f"Saved plot: {out_path}")

    plt.show()


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

    # Group by biome × taxon
    grouped = (
        df.groupby(["main_biome", taxon])[value_col]
          .mean()
          .reset_index()
        if mode == "mean"
        else df.groupby(["main_biome", taxon])[value_col]
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

    # Labels
    xlabel = LABEL_MAP.get("main_biome", "main_biome")
    ylabel = LABEL_MAP.get(value_col, value_col) if mode == "mean" else "Relative Contribution"
    title_taxon = TAXON_PLURALS.get(taxon, taxon)

    # Plot
    ax = pivot.plot(kind="bar", stacked=True, figsize=(10, 8), colormap="tab20")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if mode == "mean":
        ax.set_title(f"Metabolic Potential across Biomes (mean) by Top {top_n} {title_taxon}")
    else:
        ax.set_title(f"Relative Contribution of Top {top_n} {title_taxon} to Metabolic Potential across Biomes")

    if log_scale and mode == "mean":  # μόνο στο mean βγάζει νόημα log-scale
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

    plt.show()


# -------------------------------
# Supplementary plot 
# -------------------------------
def plot_genome_counts_by_taxon(df, taxon="phylum", top_n=20, save_dir=None, log_scale=False):
    """Plot number of genomes per phylum per biome (stacked or grouped bars)."""
    taxa = split_and_clean_taxonomy(df, "gtdb_taxonomy")
    df[taxon] = taxa[taxon]

    # Count genomes
    counts = df.groupby(["main_biome", taxon])["patric_id"].nunique().reset_index()
    top_taxa = df[taxon].value_counts().head(top_n).index
    counts = counts[counts[taxon].isin(top_taxa)]

    pivot = counts.pivot(index="main_biome", columns=taxon, values="patric_id").fillna(0)

    ax = pivot.plot(kind="bar", stacked=True, figsize=(10, 6), colormap="tab20")
    ax.set_xlabel(LABEL_MAP.get("main_biome", "Biome"))
    ax.set_ylabel("Number of genomes")

    if log_scale:
        ax.set_yscale("log")
        ax.set_title(f"Genome counts per biome for Top {top_n} {TAXON_PLURALS.get(taxon, taxon)} (log scale)")
    else:
        ax.set_title(f"Genome counts per biome for Top {top_n} {TAXON_PLURALS.get(taxon, taxon)}")

    ax.legend(
        bbox_to_anchor=(0.5, -0.25),
        loc="upper center",
        ncol=5,
        title=taxon.capitalize()
    )
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    if save_dir:
        out_path = save_dir / f"genome_counts_top{top_n}_{taxon}{'_log' if log_scale else ''}.png"
        plt.savefig(out_path, dpi=300)
        print(f"Saved plot: {out_path}")

    plt.show()


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
    metabolic_df = load_data(METABOLIC_POTENTIAL_1, filetype="csv")

    # Load metadata with biome assignment
    biomes_df = load_data(COMPACT_METADATA_WITH_BIOME)
    biomes_df = biomes_df[biomes_df["main_biome"].isin(["Soil", "Marine", "Freshwater"])]
    biomes_df = biomes_df[["patric_id", "main_biome"]]

    # Merge metabolic data with biome data
    merged = pd.merge(metabolic_df, biomes_df, on="patric_id", how="inner")

    # Save merged
    merged.to_csv("output/metabolic_potential_summary_with_biome.csv", index=False)

    # Metrics
    metrics = [
        "Total_Seeds", "Total_non_Seeds", "Ratio",
        "Seeds_per_Mbp", "Non_Seeds_per_Mbp", "Ratio_per_Mbp"
    ]

    # Output directory
    plot_dir = OUTPUT_DIR / "plots" / "biomes"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Collect stats
    kruskal_results = {}
    dunn_results = {}

    for col in metrics:
        print(f"\n=== {col} by biome ===")

        groups = [merged.loc[merged["main_biome"] == b, col].dropna()
                  for b in ["Soil", "Marine", "Freshwater"]]

        # Kruskal-Wallis
        H, p = kruskal(*groups)
        kruskal_results[col] = {"H_stat": H, "p_value": p}
        kruskal_results[col].update({b: len(g) for b, g in zip(["Soil","Marine","Freshwater"], groups)})
        print(f"Kruskal-Wallis H={H:.2f}, p={p:.2e}")

        # Dunn post-hoc
        posthoc = sp.posthoc_dunn(groups, p_adjust="bonferroni")
        posthoc.index = posthoc.columns = ["Soil", "Marine", "Freshwater"]
        dunn_results[col] = posthoc
        print(posthoc)

        # Save plot (boxplot)
        plot_boxplot_by_biome(
            merged,
            column=col,
            group_col="main_biome",
            palette="Set2",
            save_dir=plot_dir,
            log_scale=("per_Mbp" in col),
            posthoc_df=posthoc,
        )

    # ---------------------------
    # Save stats
    # ---------------------------
    save_stats(kruskal_results, plot_dir / "kruskal_results.tsv")

    for col, df in dunn_results.items():
        df.to_csv(plot_dir / f"dunn_posthoc_{col}.tsv", sep="\t")
        print(f"Saved Dunn’s post-hoc for {col}")

    # ---------------------------
    # Extra plots
    # ---------------------------
    # Stacked barplot (mean)
    plot_stacked_by_taxon(
        merged,
        value_col="Ratio_per_Mbp",
        taxon="phylum",
        top_n=15,
        save_dir=plot_dir,
        mode="mean"
    )

    # Stacked barplot (normalized)
    plot_stacked_by_taxon(
        merged,
        value_col="Ratio_per_Mbp",
        taxon="phylum",
        top_n=15,
        save_dir=plot_dir,
        mode="normalized"
    )

    # Genome counts
    plot_genome_counts_by_taxon(
        merged,
        taxon="phylum",
        top_n=15,
        log_scale=True,
        save_dir=plot_dir,
    )


if __name__ == "__main__":
    main()
