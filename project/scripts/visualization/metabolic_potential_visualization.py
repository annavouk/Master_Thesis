"""
Metabolic Potential Quantification (Approach 1) - Visualization & Statistical Analysis

Perform exploratory and statistical analysis on genome-wide metabolic potential
(measured as the ratio of seed to non-seed compounds per genome).
Visualize the data using histograms, boxplots, scatterplots and violin plots,
and assess differences between taxonomic groups with Kruskal-Wallis and Dunn's post-hoc tests.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
import scikit_posthocs as sp
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression
from scipy.stats import kruskal

from config import METABOLIC_POTENTIAL_1
from utils import (
    load_data,
    split_and_clean_taxonomy,
    plot_histogram,
    plot_boxplot,
    TAXON_PLURALS,
    LABEL_MAP,
)


# ------------------------
# Top/Bottom species by metabolic potential
# ------------------------
def top_bottom_species(df, column="Ratio_per_Mbp", top_n=10):
    """Return top and bottom species by a given metric (e.g., Ratio_per_Mbp)."""
    cols = [
        "species",
        "Ratio_per_Mbp",
        "Ratio",
        "genome_length_Mbp",
        "Total_Seeds",
        "Total_non_Seeds",
    ]

    top_df = df[cols].sort_values(column, ascending=False).head(top_n).reset_index(drop=True)
    bottom_df = df[cols].sort_values(column, ascending=True).head(top_n).reset_index(drop=True)

    print(f"\nTop {top_n} species by {column}:")
    print(top_df)
    print(f"\nBottom {top_n} species by {column}:")
    print(bottom_df)

    return top_df, bottom_df


def plot_species_metrics(df, species_col="species", save_path=None, title=None):
    """Plot Total Seeds, Total Non-Seeds for a set of species annotated with genome size."""
    species = df[species_col]
    seeds = df["Total_Seeds"]
    non_seeds = df["Total_non_Seeds"]
    genome_size = df["genome_length_Mbp"]
    ratio_per_Mbp = df["Ratio_per_Mbp"]

    x = np.arange(len(species))
    width = 0.35

    fig, ax = plt.subplots(figsize=(14, 8))

    bars1 = ax.bar(x - width/2, seeds, width, label="Total Seeds", color="teal")
    bars2 = ax.bar(x + width/2, non_seeds, width, label="Total Non-Seeds", color="lightblue")

    # Annotate bars with counts
    for i in range(len(species)):
        ax.text(x[i] - width/2, seeds[i] + max(seeds)*0.005, f"{seeds[i]}", ha="center", va="bottom", fontsize=8)
        ax.text(x[i] + width/2, non_seeds[i] + max(non_seeds)*0.005, f"{non_seeds[i]}", ha="center", va="bottom", fontsize=8)
        # Genome size annotation on top
        ax.text(x[i], max(seeds[i], non_seeds[i]) + max(seeds.max(), non_seeds.max())*0.08,
                f"{genome_size[i]:.4f} Mbp", ha="center", va="bottom", fontsize=8, color="black")
        # Ratio per Mbp annotation on top
        ax.text(x[i], max(seeds[i], non_seeds[i]) + max(seeds.max(), non_seeds.max())*0.11, 
                f"Ratio/Mbp: {ratio_per_Mbp[i]:.3f}", ha="center", va="bottom", fontsize=8, color="black")

    ax.set_xticks(x)
    ax.set_xticklabels(species, rotation=45, ha="right")
    ax.set_ylabel("Count")
    ax.set_ylim(0, max(seeds.max(), non_seeds.max()) * 1.2)
    ax.set_title(title if title else "Species Metrics", pad=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2, frameon=False)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300)
    plt.show()


# ------------------------
# Outliers
# ------------------------
def get_outliers(df, column="Ratio", taxon="phylum", top_n=10, method="2std"):
    """
    Identify genomes that are statistical outliers for a given metric,
    and summarize their taxonomic distribution.
    """
    if method == "2std":
        mean = df[column].mean()
        std = df[column].std()
        lower = mean - 2 * std
        upper = mean + 2 * std
        rule_desc = f"±2*STD (mean={mean:.3f}, std={std:.3f})"
    elif method == "iqr":
        Q1 = df[column].quantile(0.25)
        Q3 = df[column].quantile(0.75)
        IQR = Q3 - Q1
        lower = Q1 - 1.5 * IQR
        upper = Q3 + 1.5 * IQR
        rule_desc = f"IQR rule (Q1={Q1:.3f}, Q3={Q3:.3f}, IQR={IQR:.3f})"
    else:
        raise ValueError("Method must be '2std' or 'iqr'.")

    below = df[df[column] < lower]
    above = df[df[column] > upper]

    print(f"\nOutlier detection using {method.upper()} → {rule_desc}")
    print(f"Lower bound = {lower:.3f}, Upper bound = {upper:.3f}")
    print(f"Genomes BELOW: {len(below)} / {len(df)} ({100*len(below)/len(df):.2f}%)")
    print(f"Genomes ABOVE: {len(above)} / {len(df)} ({100*len(above)/len(df):.2f}%)")

    print(f"\nTop {top_n} {taxon} BELOW bound:")
    print(below[taxon].value_counts().head(top_n))
    print(f"\nTop {top_n} {taxon} ABOVE bound:")
    print(above[taxon].value_counts().head(top_n))

    return {
        "below": below,
        "above": above,
        "lower": lower,
        "upper": upper,
        "taxon": taxon,
        "column": column,
        "method": method,
        "rule": rule_desc,
    }


def plot_outliers(
    df,
    outlier_res,
    top_n=10,
    color_below="blue",
    color_above="red",
    title=None,
    xlabel=None,
    ylabel=None,
    save_plot=False,
    plot_path=None,
):
    """Barplots of the distribution of outlier genomes by taxonomic group."""
    below = outlier_res["below"]
    above = outlier_res["above"]
    taxon = outlier_res["taxon"]
    column = outlier_res["column"]
    method = outlier_res["method"]
    rule_desc = outlier_res["rule"]

    plural_taxon = TAXON_PLURALS.get(taxon, taxon)

    below_counts = below[taxon].value_counts().head(top_n)
    above_counts = above[taxon].value_counts().head(top_n)
    total_counts = df[taxon].value_counts()

    fig, axs = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    default_title = f"Distribution of Outlier Genomes by {column}"
    fig.suptitle(title if title is not None else default_title, fontsize=15, y=0.98)
    fig.supxlabel(
        xlabel if xlabel is not None else f"Top {top_n} {plural_taxon.capitalize()}",
        fontsize=12,
        y=0.1,
    )

    # Below plot
    bars = axs[0].bar(
        below_counts.index, below_counts.values, color=color_below, edgecolor="black"
    )
    axs[0].set_ylabel(ylabel if ylabel else "Number of Genomes")
    axs[0].tick_params(axis="x", rotation=45)
    axs[0].text(
        0.98,
        0.92,
        f"N = {len(below):,} genomes",
        transform=axs[0].transAxes,
        ha="right",
        va="top",
        color=color_below,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=color_below),
    )
    for i, bar in enumerate(bars):
        height = bar.get_height()
        taxon_name = below_counts.index[i]
        total = total_counts[taxon_name]
        axs[0].annotate(
            f"{int(height)}/{total}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7,
        )

    # Above plot
    bars = axs[1].bar(
        above_counts.index, above_counts.values, color=color_above, edgecolor="black"
    )
    axs[1].tick_params(axis="x", rotation=45)
    axs[1].text(
        0.98,
        0.92,
        f"N = {len(above):,} genomes",
        transform=axs[1].transAxes,
        ha="right",
        va="top",
        color=color_above,
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=color_above),
    )
    for i, bar in enumerate(bars):
        height = bar.get_height()
        taxon_name = above_counts.index[i]
        total = total_counts[taxon_name]
        axs[1].annotate(
            f"{int(height)}/{total}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7,
        )

    # Legend dynamic based on method
    legend_elems = [
        Patch(facecolor=color_below, edgecolor="black", label=f"Below {rule_desc}"),
        Patch(facecolor=color_above, edgecolor="black", label=f"Above {rule_desc}"),
    ]
    fig.legend(
        handles=legend_elems, loc="lower center", ncol=2, frameon=False, fontsize=12
    )

    if save_plot and plot_path:
        fig.savefig(plot_path, dpi=300)
    fig.tight_layout()
    plt.show()
    return fig


# ------------------------
# Regression
# ------------------------
def plot_regression_line(ax, xvals, yvals, color="red"):
    """Fit and plot a regression line on the given axes."""
    model = LinearRegression().fit(xvals.reshape(-1, 1), yvals)
    x_range = np.linspace(xvals.min(), xvals.max(), 100)
    y_pred = model.predict(x_range.reshape(-1, 1))
    ax.plot(x_range, y_pred, color=color, linewidth=2, label="Regression line")


def correlation_stats(
    df, x="genome_length_Mbp", y="Ratio_per_Mbp", log_x=False, log_y=True
):
    """Calculate and print Pearson and Spearman correlation coefficients"""
    vals_x = df[x].dropna()
    vals_y = df[y].dropna()

    if log_x:
        vals_x = np.log10(vals_x[vals_x > 0])
    if log_y:
        vals_y = np.log10(vals_y[vals_y > 0])

    ix = vals_x.index.intersection(vals_y.index)
    vals_x = vals_x.loc[ix]
    vals_y = vals_y.loc[ix]
    # Pearson
    pearson, p_pearson = pearsonr(vals_x, vals_y)
    # Spearman
    spearman, p_spearman = spearmanr(vals_x, vals_y)
    print(f"Pearson r = {pearson:.3f} (p={p_pearson:.2e})")
    print(f"Spearman r = {spearman:.3f} (p={p_spearman:.2e})")
    return pearson, p_pearson, spearman, p_spearman


# ------------------------
# Scatter plots
# ------------------------
def plot_scatter(
    df,
    x="genome_length_Mbp",
    y="Ratio_per_Mbp",
    hue=None,
    taxon="phylum",
    top_n=10,
    title=None,
    xlabel=None,
    ylabel=None,
    legend_title=None,
    palette=None,
    xscale="linear",
    yscale="log",
    save_plot=False,
    plot_path=None,
):
    """
    Scatter plot for two metrics with optional log/linear scaling, regression line and
    correlation statistics. Highlights top taxa if specified.
    """
    top_taxa = df[taxon].value_counts().nlargest(top_n).index
    sub = df[df[taxon].isin(top_taxa)]
    filtered = sub[[x, y, taxon]].dropna()

    if xscale != "linear":
        filtered = filtered[filtered[x] > 0]
    if yscale != "linear":
        filtered = filtered[filtered[y] > 0]

    n = len(filtered)
    taxon_plural = TAXON_PLURALS.get(taxon, taxon)

    title = (
        title
        or f"Variation of {LABEL_MAP.get(y, y)} with {LABEL_MAP.get(x, x)} among top {top_n} {taxon_plural}\n(N = {n} genomes)"
    )
    xlabel = xlabel or LABEL_MAP.get(x, x)
    ylabel = ylabel or LABEL_MAP.get(y, y)
    legend_title = legend_title or taxon_plural.capitalize()

    if xscale != "linear":
        xlabel += " (log scale)"
    if yscale != "linear":
        ylabel += " (log scale)"

    plot_hue = hue if hue is not None else taxon

    fig, ax = plt.subplots(figsize=(15, 6))
    sns.scatterplot(
        data=filtered, x=x, y=y, hue=plot_hue, alpha=0.7, ax=ax, palette=palette
    )

    # Prepare regression (use log-transform if needed)
    xv = filtered[x].values
    yv = filtered[y].values
    logx = xscale != "linear"
    logy = yscale != "linear"
    plot_xv = np.log10(xv) if logx else xv
    plot_yv = np.log10(yv) if logy else yv

    # Fit regression line and plot
    if len(plot_xv) > 1 and len(plot_yv) > 1:
        model = LinearRegression().fit(plot_xv.reshape(-1, 1), plot_yv)
        x_range = np.linspace(plot_xv.min(), plot_xv.max(), 100)
        y_pred = model.predict(x_range.reshape(-1, 1))
        # For display, invert the log transform if axis is log
        if logy:
            y_plot = 10**y_pred
        else:
            y_plot = y_pred
        if logx:
            x_plot = 10**x_range
        else:
            x_plot = x_range
        ax.plot(x_plot, y_plot, color="red", linewidth=2, label="Regression line")

        # Correlations
        pearson, p1 = pearsonr(plot_xv, plot_yv)
        spearman, p2 = spearmanr(plot_xv, plot_yv)

        threshold = 1e-16
        if p1 < threshold:
            p1_str = f"p < {threshold:.0e}"
        else:
            p1_str = f"p = {p1:.2e}"
        if p2 < threshold:
            p2_str = f"p < {threshold:.0e}"
        else:
            p2_str = f"p = {p2:.2e}"

        # Annotation
        ax.annotate(
            f"Pearson r = {pearson:.2f}, {p1_str}\nSpearman r = {spearman:.2f}, {p2_str}",
            xy=(0.98, 0.98),
            xycoords="axes fraction",
            ha="right",
            va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9),
            fontsize=11,
        )

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", title=legend_title)
    if save_plot and plot_path:
        fig.savefig(plot_path, dpi=300)
    fig.tight_layout()
    plt.show()
    return fig


def plot_scatter_continuous_colormap(
    df,
    x,
    y,
    color_col=None,
    title=None,
    xlabel=None,
    ylabel=None,
    xscale="linear",
    yscale="linear",
    cmap="viridis",
    regression=True,
    save_plot=False,
    plot_path=None,
):
    """Plot a scatter plot with continuous color mapping for a third variable."""
    filtered = df[[x, y, color_col] if color_col else [x, y]].dropna()
    xv = filtered[x].values
    yv = filtered[y].values
    cvals = filtered[color_col].values if color_col else None

    fig, ax = plt.subplots(figsize=(12, 6))
    sc = (
        ax.scatter(xv, yv, c=cvals, cmap=cmap, alpha=0.6, s=15)
        if color_col
        else ax.scatter(xv, yv, alpha=0.7, s=15)
    )
    if color_col:
        plt.colorbar(sc, ax=ax, label=color_col.replace("_", " ").capitalize())

    if regression and len(xv) > 1 and len(yv) > 1:
        model = LinearRegression().fit(xv.reshape(-1, 1), yv)
        x_range = np.linspace(xv.min(), xv.max(), 100).reshape(-1, 1)
        y_pred = model.predict(x_range)
        ax.plot(x_range, y_pred, color="red", lw=2, label="Regression line")
        pearson, p1 = pearsonr(xv, yv)
        spearman, p2 = spearmanr(xv, yv)

        threshold = 1e-16
        if p1 < threshold:
            p1_str = f"p < {threshold:.0e}"
        else:
            p1_str = f"p = {p1:.2e}"
        if p2 < threshold:
            p2_str = f"p < {threshold:.0e}"
        else:
            p2_str = f"p = {p2:.2e}"

        ax.annotate(
            f"Pearson r = {pearson:.2f}, {p1_str}\nSpearman r = {spearman:.2f}, {p2_str}",
            xy=(0.05, 0.95),
            xycoords="axes fraction",
            ha="left",
            va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9),
            fontsize=11,
        )

    ax.set_title(title or f"{y} vs {x}")
    ax.set_xlabel(xlabel or x)
    ax.set_ylabel(ylabel or y)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if save_plot and plot_path:
        fig.savefig(plot_path, dpi=300)
    fig.tight_layout()
    plt.show()
    return fig


# ------------------------
# Violin plots
# ------------------------
def plot_violin_ratio_by_taxon(
    df,
    taxon="phylum",
    column="Ratio_per_Mbp",
    top_n=10,
    y_limit=(0, 5),
    save_plot=False,
    plot_path=None,
):
    """Generate violin plots for total number of seeds to non-seeds per phylum."""
    top_taxa = df[taxon].value_counts().nlargest(top_n).index
    df_sub = df[df[taxon].isin(top_taxa)]

    n_per_taxon = df_sub[taxon].value_counts().reindex(top_taxa)
    avg_size_per_taxon = (
        df_sub.groupby(taxon)["genome_length_Mbp"]
        .mean()
        .reindex(top_taxa)
        .round(0)
        .astype(int)
    )

    medians = df_sub.groupby(taxon)[column].median().reindex(top_taxa)

    label = LABEL_MAP.get(column, column)

    sns.set(style="whitegrid")

    fig, ax = plt.subplots(figsize=(14, 6))
    sns.violinplot(
        data=df_sub,
        x=taxon,
        y=column,
        inner="box",
        palette="pastel",
        legend=False,
        order=top_taxa,
        ax=ax,
    )
    plt.title(f"Distribution of {label} per {taxon.capitalize()}", fontsize=16)
    plt.xticks(rotation=45, ha="right")

    xticklabels = [
        f"{taxon}\n(N={n_per_taxon[taxon]:,}, avg size={avg_size_per_taxon[taxon]:,} Mbp)"
        for taxon in top_taxa
    ]
    ax.set_xticklabels(xticklabels)

    for i, phylum in enumerate(top_taxa):
        ax.scatter(
            i,
            medians[phylum],
            color="black",
            marker="o",
            s=50,
            zorder=10,
            label="Median" if i == 0 else "",
        )

    # Set the y-axis limit
    plt.ylim(y_limit)
    if save_plot and plot_path:
        fig.savefig(plot_path, dpi=300)
    fig.tight_layout()
    plt.show()
    return fig


# ------------------------
# Stats test
# ------------------------
def kruskal_test_by_taxon(df, column, taxon, top_n=10):
    """Perform Kruskal-Wallis test for a given column by taxonomic group."""
    top_taxa = df[taxon].value_counts().nlargest(top_n).index
    sub = df[df[taxon].isin(top_taxa)]
    data = [sub[sub[taxon] == group][column].dropna() for group in top_taxa]
    H, p = kruskal(*data)

    threshold = 1e-16
    if p < threshold:
        print(f"Kruskal–Wallis for {column} by {taxon}: H={H:.2f}, p < {threshold:.0e}")
    else:
        print(f"Kruskal–Wallis for {column} by {taxon}: H={H:.2f}, p = {p:.2e}")

    return H, p


def dunn_posthoc_test(
    df,
    taxon="phylum",
    column="Ratio_per_Mbp",
    top_n=10,
    p_adjust="fdr_bh",
    export_csv=None,
    alpha=0.05,
    verbose=True,
):
    """Perform Dunn's post-hoc test after Kruskal-Wallis for a given column by taxonomic group."""
    if verbose:
        # Subset only the top_n groups
        top_taxa = df[taxon].value_counts().nlargest(top_n).index
        sub = df[df[taxon].isin(top_taxa)]

    # Run Dunn’s post-hoc test
    dunn_df = sp.posthoc_dunn(sub, val_col=column, group_col=taxon, p_adjust=p_adjust)

    if export_csv:
        dunn_df.to_csv(export_csv)
        if verbose:
            print(f"[INFO] Dunn’s test results exported to {export_csv}")

    # Extract significant results
    significant = dunn_df < alpha
    sig_pairs = []
    for i in significant.index:
        for j in significant.columns:
            if i != j and significant.loc[i, j]:
                pair = tuple(sorted((i, j)))
                sig_pairs.append((pair[0], pair[1], dunn_df.loc[i, j]))
    # Remove duplicate pairs
    sig_pairs = list({(a, b): p for (a, b, p) in sig_pairs}.items())

    if verbose:
        print(f"Significant comparisons (p < {alpha}): {len(sig_pairs)} unique pairs.")
        for (a, b), p in sig_pairs:
            print(f"{a} vs {b}: p = {p:.2e}")

    return dunn_df


# ------------------------
# Main
# ------------------------
def main():
    # Load and split taxonomy
    df = load_data(METABOLIC_POTENTIAL_1, filetype="csv")
    taxonomy_df = split_and_clean_taxonomy(df, "gtdb_taxonomy")
    df = pd.concat([df, taxonomy_df], axis=1)

    # Species with the highest and lowest metabolic potential
    top_df, bottom_df = top_bottom_species(df, column="Ratio_per_Mbp", top_n=10)
    plot_species_metrics(
        top_df,
        species_col="species",
        save_path="output/plots/metabolic_potential_1/species_metrics_top_ratio_Mbp.png",
        title="Metabolic Potential Metrics Across Top Species by Seed/Non-Seed per Genome Size Ratio",
    )
    plot_species_metrics(
        bottom_df,
        species_col="species",
        save_path="output/plots/metabolic_potential_1/species_metrics_bottom_ratio_Mbp.png",
        title="Metabolic Potential Metrics Across Bottom Species by Seed/Non-Seed per Genome Size Ratio",
    )

    stats_file = "output/metabolic_potential_summary_stats.txt"

    # Clear file for summary statistics
    open(stats_file, "w").close()

    # Plot and save distributions (raw and normalized) for Seeds, non-Seeds, and their ratio.
    plot_histogram(
        df,
        column="Total_Seeds",
        title="Distribution of Total Seeds per Genome",
        xlabel="Total number of Seeds per Genome",
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/histogram_total_seeds.png",
        save_stats=True,
        stats_path=stats_file,
    )
    plot_histogram(
        df,
        column="Total_non_Seeds",
        title="Distribution of Total Non-Seeds per Genome",
        xlabel="Total number of non-Seeds per Genome",
        color="skyblue",
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/histogram_total_non_seeds.png",
        save_stats=True,
        stats_path=stats_file,
    )
    plot_histogram(
        df,
        column="Ratio",
        title="Distribution of Total Seeds/Non-Seeds Ratio per Genome",
        xlabel="Seeds/non-Seeds Ratio per Genome",
        color="orange",
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/histogram_ratio.png",
        save_stats=True,
        stats_path=stats_file,
    )
    plot_histogram(
        df,
        column="Seeds_per_Mbp",
        title="Distribution of Total Seeds per Genome size per Genome",
        xlabel="Total number of Seeds/Genome Size (Mbp) per Genome (log scale)",
        color="green",
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/histogram_total_seeds_per_Mbp.png",
        save_stats=True,
        stats_path=stats_file,
        log_scale=True,
    )
    plot_histogram(
        df,
        column="Non_Seeds_per_Mbp",
        title="Distribution of Total Non-Seeds per Genome size per Genome",
        xlabel="Total number of non-Seeds/Genome Size (Mbp) per Genome (log scale)",
        color="lime",
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/histogram_total_non_seeds_per_Mbp.png",
        save_stats=True,
        stats_path=stats_file,
        log_scale=True,
    )
    plot_histogram(
        df,
        column="Ratio_per_Mbp",
        title="Distribution of Seeds/Non-Seeds Ratio per Genome size per Genome",
        xlabel="Seeds/non-Seeds/Genome Size (Mbp) per Genome (log scale)",
        color="purple",
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/histogram_ratio_per_Mbp.png",
        save_stats=True,
        stats_path=stats_file,
        log_scale=True,
    )

    # Plot boxplots for seeds/non-seeds and seeds/non-seeds/genome_size ratios.
    plot_boxplot(
        df,
        column="Ratio",
        title="Seeds/Non-Seeds Ratio per Genome",
        xlabel="Seeds/Non-Seeds Ratio",
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/boxplot_ratio.png",
        method="2std",    
    )
    plot_boxplot(
        df,
        column="Ratio_per_Mbp",
        color="lightblue",
        title="Seeds/Non-Seeds Ratio per Genome Size per Genome",
        xlabel="Seeds/Non-Seeds Ratio/Genome Size (Mbp) per Genome (log scale)",
        log_scale=True,
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/boxplot_ratio_per_Mbp.png",
        method="iqr",
    )
    
    # Identify and plot outliers for seeds/non-seeds (raw and normalized by genome size) by phylum.
    outliers_ratio = get_outliers(df, column="Ratio", taxon="phylum", top_n=10, method="2std")
    plot_outliers(df, outliers_ratio,
              save_plot=True,
              plot_path="output/plots/metabolic_potential_1/outliers_phylum_ratio.png")

    outliers_mbp = get_outliers(df, column="Ratio_per_Mbp", taxon="phylum", top_n=10, method="iqr")
    plot_outliers(df, outliers_mbp,
              save_plot=True,
              plot_path="output/plots/metabolic_potential_1/outliers_phylum_ratio_per_Mbp.png")

    # Plot scatter plot for seed/non-seed/genome_size vs genome size
    plot_scatter(
        df,
        title="Variation of Seeds to Non-Seeds to Genome Size with Genome Size among Top 10 Phyla",
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/scatter_ratio_per_Mbp_vs_genome_size.png",
    )
    # Plot scatter plot for seed/non-seed vs genome size
    plot_scatter(
        df,
        title="Variation of Seeds to Non-Seeds with Genome Size among Top 10 Phyla",
        y="Ratio",
        yscale="linear",
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/scatter_ratio_vs_genome_size.png",
    )

    # Plot scatter plot for seed vs non-seed
    plot_scatter_continuous_colormap(
        df,
        "Total_Seeds",
        "Total_non_Seeds",
        color_col="genome_length_Mbp",
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/scatter_seed_vs_non_seed.png",
    )

    # Plot violin plot for seeds/non-seeds ratio/genome size (Mbp) by phylum
    plot_violin_ratio_by_taxon(
        df,
        top_n=10,
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/violin_ratio_Mbp_phyla.png",
    )

    # Plot violin plot for seeds/non-seeds ratio by phylum
    plot_violin_ratio_by_taxon(
        df,
        column="Ratio",
        top_n=10,
        y_limit=(0.2, 0.6),
        save_plot=True,
        plot_path="output/plots/metabolic_potential_1/violin_ratio_phyla.png",
    )

    # Kruskal-Wallis test for Ratio by phylum
    print("\nKruskal–Wallis test for Ratio per phylum:")
    kruskal_test_by_taxon(df, "Ratio_per_Mbp", "phylum")

    # Dunn's post-hoc for Ratio by phylum
    dunn_df = dunn_posthoc_test(df, taxon="phylum", column="Ratio_per_Mbp", top_n=10)
    print("\nDunn's post-hoc test results:")
    print(dunn_df.head())


if __name__ == "__main__":
    main()
