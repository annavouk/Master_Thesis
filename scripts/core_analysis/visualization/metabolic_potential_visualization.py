"""
Metabolic Potential Quantification - Visualization & Statistical Analysis

Perform exploratory and statistical analysis on genome-wide metabolic potential
(measured as the ratio of seed to non-seed compounds per genome).
Visualize the data using histograms, scatterplots, violin plots and barplots
and assess differences between taxonomic groups with Kruskal-Wallis and Dunn's post-hoc tests.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
import scikit_posthocs as sp

from scipy.stats import pearsonr, spearmanr, kruskal
from sklearn.linear_model import LinearRegression

from config import (
    METABOLIC_POTENTIAL_TSV,  # input
    METABOLIC_POTENTIAL_PLOTS_DIR,  # plots dir
    OUTPUT_DIR,  # output dir
)
from utils import (
    load_data,
    split_and_clean_taxonomy,
    plot_histogram,
    TAXON_PLURALS,
    LABEL_MAP,
)


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

    print(f"\nOutlier detection using {method.upper()} -> {rule_desc}")
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
    rule_desc = outlier_res["rule"]

    plural_taxon = TAXON_PLURALS.get(taxon, taxon)

    if outlier_res["below"].empty and outlier_res["above"].empty:
        print("No outliers detected; skipping plot_outliers().")
        return None

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
    fig.subplots_adjust(bottom=0.3)

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
        total = total_counts.get(taxon_name, 0)
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
        total = total_counts.get(taxon_name, 0)
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
        fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.tight_layout()
    plt.close(fig)
    return None


# ------------------------
# Top/Bottom species by metabolic potential
# ------------------------
def top_bottom_species(df, column="Ratio", top_n=10):
    """Return species with the highest/lowest values of a given metric (e.g. Ratio)."""
    cols = [
        "species",
        "Ratio",
        "genome_size_Mbp",
        "Total_Seeds",
        "Total_non_Seeds",
        "Seeds_per_Mbp",
        "Non_Seeds_per_Mbp",
    ]

    top_df = (
        df[cols].sort_values(column, ascending=False).head(top_n).reset_index(drop=True)
    )
    bottom_df = (
        df[cols].sort_values(column, ascending=True).head(top_n).reset_index(drop=True)
    )

    print(f"\nTop {top_n} species by {column}:")
    print(top_df)
    print(f"\nBottom {top_n} species by {column}:")
    print(bottom_df)

    return top_df, bottom_df


def plot_species_metrics(
    df,
    species_col="species",
    save_path=None,
    title=None,
    ratio_col="Ratio",
    seeds_col="Total_Seeds",
    non_seeds_col="Total_non_Seeds",
    genome_size_col="genome_size_Mbp",
):
    """Plot Total Seeds, Total Non-Seeds for a set of species annotated with genome size and ratio."""
    species = df[species_col].astype(str).to_numpy()
    seeds = pd.to_numeric(df[seeds_col], errors="coerce").to_numpy()
    non_seeds = pd.to_numeric(df[non_seeds_col], errors="coerce").to_numpy()

    ratio_vals = (
        pd.to_numeric(df[ratio_col], errors="coerce").to_numpy()
        if ratio_col in df.columns
        else None
    )
    genome_size = (
        pd.to_numeric(df[genome_size_col], errors="coerce").to_numpy()
        if genome_size_col and genome_size_col in df.columns
        else None
    )

    seeds = np.nan_to_num(seeds, nan=0.0)
    non_seeds = np.nan_to_num(non_seeds, nan=0.0)

    x = np.arange(len(species))
    width = 0.35
    fig, ax = plt.subplots(figsize=(14, 8))

    ax.bar(x - width / 2, seeds, width, label="Total Seeds", color="teal")
    ax.bar(x + width / 2, non_seeds, width, label="Total Non-Seeds", color="lightblue")

    max_bar = float(max(seeds.max(), non_seeds.max(), 1.0))
    ax.set_ylim(0, max_bar * 1.22)
    small_off = max_bar * 0.01
    gsize_off = max_bar * 0.06
    ratio_off = max_bar * 0.09

    for i in range(len(species)):
        ax.text(
            x[i] - width / 2,
            seeds[i] + small_off,
            f"{int(seeds[i])}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
        ax.text(
            x[i] + width / 2,
            non_seeds[i] + small_off,
            f"{int(non_seeds[i])}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

        top_here = max(seeds[i], non_seeds[i])
        if genome_size is not None and np.isfinite(genome_size[i]):
            ax.text(
                x[i],
                top_here + gsize_off,
                f"{genome_size[i]:.4f} Mbp",
                ha="center",
                va="bottom",
                fontsize=8,
                color="black",
            )
        if ratio_vals is not None and np.isfinite(ratio_vals[i]):
            ax.text(
                x[i],
                top_here + ratio_off,
                f"Ratio: {ratio_vals[i]:.3f}",
                ha="center",
                va="bottom",
                fontsize=8,
                color="black",
            )

    ax.set_xticks(x)
    ax.set_xticklabels(species, rotation=45, ha="right")
    ax.set_ylabel("Count")
    ax.set_title(title if title else "Species Metrics", pad=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.4), ncol=2, frameon=False)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return None


# ------------------------
# Violin plots
# ------------------------
def plot_violin_ratio_by_taxon(
    df,
    taxon="phylum",
    column="Ratio",
    top_n=10,
    y_limit=(0, 5),
    save_plot=False,
    plot_path=None,
):
    """Generate violin plots of a specific column (e.g. Ratio) across the top-N groups of a taxonomic rank (e.g. phylum)."""
    top_taxa = df[taxon].value_counts().nlargest(top_n).index
    df_sub = df[df[taxon].isin(top_taxa)]

    n_per_taxon = df_sub[taxon].value_counts().reindex(top_taxa)

    avg_size_per_taxon = (
        df_sub.groupby(taxon)["genome_size_Mbp"]
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
        order=top_taxa,
        ax=ax,
    )
    plt.title(f"Distribution of {label} per {taxon.capitalize()}", fontsize=16)
    plt.xticks(rotation=45, ha="right")
    plt.subplots_adjust(bottom=0.5)

    xticklabels = [
        f"{tx}\n(N={n_per_taxon[tx]:,}, avg size={avg_size_per_taxon[tx]:,} Mbp)"
        for tx in top_taxa
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
        fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.tight_layout()
    plt.close(fig)
    return None


# ------------------------
# Scatter plots
# ------------------------
def plot_scatter(
    df,
    x="genome_size_Mbp",
    y="Ratio",
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

    fig, ax = plt.subplots(figsize=(14, 6))
    sns.scatterplot(
        data=filtered, x=x, y=y, hue=plot_hue, alpha=0.7, ax=ax, palette=palette
    )

    # Prepare regression (use log-transform if needed)
    if n >= 3:
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
            x_plot = 10**x_range if logx else x_range
            y_plot = 10**y_pred if logy else y_pred
            ax.plot(x_plot, y_plot, color="red", linewidth=2, label="Regression line")

            pearson_r, p1 = pearsonr(plot_xv, plot_yv)
            spearman_r, p2 = spearmanr(plot_xv, plot_yv)
            threshold = 1e-16
            p1_str = f"p < {threshold:.0e}" if p1 < threshold else f"p = {p1:.2e}"
            p2_str = f"p < {threshold:.0e}" if p2 < threshold else f"p = {p2:.2e}"

            ax.annotate(
                f"Pearson r = {pearson_r:.2f}, {p1_str}\nSpearman r = {spearman_r:.2f}, {p2_str}",
                xy=(0.98, 0.98),
                xycoords="axes fraction",
                ha="right",
                va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9),
                fontsize=11,
            )
        else:
            print("[plot_scatter] Not enough variation for regression/annotation.")
    else:
        print(
            f"[plot_scatter] Too few points after filtering (N={n}). Skipping regression/annotation."
        )

    # Axis/legend/save
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left", title=legend_title, fontsize=7)
    if save_plot and plot_path:
        fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.tight_layout()
    plt.close(fig)
    return None


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
    # Subframe + NaN drop
    cols = [x, y] + ([color_col] if color_col else [])
    filtered = df[cols].dropna()

    # Log-scale filter
    if xscale != "linear":
        filtered = filtered[filtered[x] > 0]
    if yscale != "linear":
        filtered = filtered[filtered[y] > 0]

    if filtered.empty:
        print(
            "[plot_scatter_continuous_colormap] No data after filtering; skipping plot."
        )
        return None

    xv = filtered[x].to_numpy()
    yv = filtered[y].to_numpy()
    cvals = filtered[color_col].to_numpy() if color_col else None

    # Labels/titles
    lx = LABEL_MAP.get(x, x)
    ly = LABEL_MAP.get(y, y)
    lc = LABEL_MAP.get(color_col, color_col) if color_col else None
    title = title or f"{ly} vs {lx}"
    xlabel = (xlabel or lx) + (" (log scale)" if xscale != "linear" else "")
    ylabel = (ylabel or ly) + (" (log scale)" if yscale != "linear" else "")

    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    if color_col:
        sc = ax.scatter(xv, yv, c=cvals, cmap=cmap, alpha=0.6, s=15)
        clabel = (lc or color_col).replace("_", " ").capitalize()
        plt.colorbar(sc, ax=ax, label=clabel)
    else:
        ax.scatter(xv, yv, alpha=0.7, s=15)

    # Regression
    if regression and xv.size >= 3 and yv.size >= 3:
        # Fit by scales
        rx = np.log10(xv) if xscale != "linear" else xv
        ry = np.log10(yv) if yscale != "linear" else yv

        if np.unique(rx).size > 1 and np.unique(ry).size > 1:
            model = LinearRegression().fit(rx.reshape(-1, 1), ry)
            xr = np.linspace(rx.min(), rx.max(), 100)
            yp = model.predict(xr.reshape(-1, 1))

            x_plot = 10**xr if xscale != "linear" else xr
            y_plot = 10**yp if yscale != "linear" else yp
            ax.plot(x_plot, y_plot, color="red", lw=2, label="Regression line")

            # Correlation stats
            pr, p1 = pearsonr(rx, ry)
            sr, p2 = spearmanr(rx, ry)
            thr = 1e-16
            p1s = f"p < {thr:.0e}" if p1 < thr else f"p = {p1:.2e}"
            p2s = f"p < {thr:.0e}" if p2 < thr else f"p = {p2:.2e}"
            ax.annotate(
                f"Pearson r = {pr:.2f}, {p1s}\nSpearman r = {sr:.2f}, {p2s}",
                xy=(0.05, 0.95),
                xycoords="axes fraction",
                ha="left",
                va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9),
                fontsize=11,
            )
        else:
            print(
                "[plot_scatter_continuous_colormap] Not enough variation for regression/annotation."
            )
    elif regression:
        print(
            f"[plot_scatter_continuous_colormap] Too few points (N={len(filtered)}). Skipping regression/annotation."
        )

    # Axes, scales, save
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    if save_plot and plot_path:
        fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.tight_layout()
    plt.close(fig)
    return None


# ------------------------
# Heatmap
# ------------------------
def plot_dunn_heatmap(dunn_df, save_path=None, figsize=(10, 8), annot=False, fmt=".1f"):
    """Plot a heatmap of Dunn's post-hoc test results using -log10(p-values)."""
    log_p = -np.log10(dunn_df.replace(0, 1e-300))
    mask = np.tril(np.ones_like(log_p, dtype=bool))
    fig = plt.figure(figsize=figsize)
    sns.heatmap(
        log_p,
        mask=mask,
        cmap="viridis",
        annot=annot,  # show actual -log10(p) values if annot=True
        fmt=fmt if annot else "",  # format for annotations
        cbar_kws={"label": r"-log$_{10}$(p-value)"},
        xticklabels=True,
        yticklabels=True,
    )
    plt.title("Dunn's Post-hoc Test (-log10 p-values)", fontsize=14)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return None


# ------------------------
# Correlation stats
# ------------------------
def correlation_stats(
    df, x="genome_size_Mbp", y="Ratio", log_x=False, log_y=False, p_threshold=1e-16
):
    """Calculate and print Pearson and Spearman correlation between two variables."""
    mask = df[[x, y]].notna().all(axis=1)
    if log_x:
        mask &= df[x] > 0
    if log_y:
        mask &= df[y] > 0

    vals_x = df.loc[mask, x].to_numpy(dtype=float)
    vals_y = df.loc[mask, y].to_numpy(dtype=float)

    if log_x:
        vals_x = np.log10(vals_x)
    if log_y:
        vals_y = np.log10(vals_y)

    if vals_x.size < 3 or np.unique(vals_x).size < 2 or np.unique(vals_y).size < 2:
        print(
            "Pearson/Spearman: not enough variation or too few points after filtering."
        )
        return np.nan, np.nan, np.nan, np.nan

    pearson_r, p_pearson = pearsonr(vals_x, vals_y)
    spearman_r, p_spearman = spearmanr(vals_x, vals_y)

    p1 = f"< {p_threshold:.0e}" if p_pearson < p_threshold else f"= {p_pearson:.2e}"
    p2 = f"< {p_threshold:.0e}" if p_spearman < p_threshold else f"= {p_spearman:.2e}"
    print(f"Pearson r = {pearson_r:.3f} (p {p1})")
    print(f"Spearman r = {spearman_r:.3f} (p {p2})")
    return pearson_r, p_pearson, spearman_r, p_spearman


def kruskal_test_by_taxon(df, column, taxon, top_n=10):
    """Perform Kruskal-Wallis test for a given column by taxonomic group."""
    top_taxa = df[taxon].value_counts().nlargest(top_n).index
    sub = df[df[taxon].isin(top_taxa)]

    groups, data = [], []
    for g in top_taxa:
        arr = sub.loc[sub[taxon] == g, column].dropna().values
        if arr.size > 0:
            groups.append(g)
            data.append(arr)

    if len(data) < 2:
        print(f"Kruskal–Wallis for {column} by {taxon}: not enough non-empty groups.")
        return np.nan, np.nan

    H, p = kruskal(*data)  # unpacks the list
    threshold = 1e-16
    msg_p = f"p < {threshold:.0e}" if p < threshold else f"p = {p:.2e}"
    print(f"Kruskal–Wallis for {column} by {taxon}: H={H:.2f}, {msg_p} (k={len(data)})")
    return H, p


def dunn_posthoc_test(
    df,
    taxon="phylum",
    column="Ratio",
    top_n=10,
    p_adjust="fdr_bh",
    export_tsv=None,
    alpha=0.05,
    verbose=True,
):
    """Perform Dunn's post-hoc test for a given column by taxonomic group."""
    # Subset only the top_n groups
    top_taxa = df[taxon].value_counts().nlargest(top_n).index
    sub = df[df[taxon].isin(top_taxa)]

    # Keep rows with non-missing target values
    sub = sub.loc[sub[column].notna()]
    # Must have at least two groups with data
    groups_present = sub[taxon].value_counts()
    if (groups_present > 0).sum() < 2:
        if verbose:
            print(
                f"Dunn post-hoc for {column} by {taxon}: not enough groups with data."
            )
        return pd.DataFrame()

    # Run Dunn’s post-hoc test (via scikit-posthocs)
    dunn_df = sp.posthoc_dunn(sub, val_col=column, group_col=taxon, p_adjust=p_adjust)

    if export_tsv:
        dunn_df.to_csv(export_tsv, sep="\t")
        if verbose:
            print(f"[INFO] Dunn’s test results exported to {export_tsv}")

    # Identify significant pairs once (use upper triangle to avoid duplicates)
    sig_mask = dunn_df.values < alpha
    upper = np.triu(np.ones_like(sig_mask, dtype=bool), k=1)
    sig_upper = sig_mask & upper

    if verbose:
        i_idx, j_idx = np.where(sig_upper)
        pairs = []
        for i, j in zip(i_idx, j_idx):
            g1 = dunn_df.index[i]
            g2 = dunn_df.columns[j]
            pval = dunn_df.iloc[i, j]
            pairs.append((str(g1), str(g2), float(pval)))

        print(f"Significant comparisons (p < {alpha}): {len(pairs)} unique pairs.")
        for a, b, pval in pairs:
            print(f"{a} vs {b}: p = {pval:.2e}")

    return dunn_df


# ------------------------
# Main
# ------------------------
def main():
    # Ensure output directories exist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    METABOLIC_POTENTIAL_PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    # Load and split taxonomy
    df = load_data(METABOLIC_POTENTIAL_TSV, filetype="tsv")
    taxonomy_df = split_and_clean_taxonomy(df, "gtdb_taxonomy")
    df = pd.concat([df, taxonomy_df], axis=1)

    stats_file = OUTPUT_DIR / "metabolic_potential_summary_stats.txt"
    
    # Clear file for summary statistics
    open(stats_file, "w").close()

    # Plot and save distributions (raw and normalized) for Seeds, non-Seeds, and their ratio.
    plot_histogram(
        df,
        column="Total_Seeds",
        title="Distribution of Total Seeds per Genome",
        xlabel="Total number of Seeds per Genome",
        save_plot=True,
        plot_path=METABOLIC_POTENTIAL_PLOTS_DIR / "histogram_total_seeds.png",
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
        plot_path=METABOLIC_POTENTIAL_PLOTS_DIR / "histogram_total_non_seeds.png",
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
        plot_path=METABOLIC_POTENTIAL_PLOTS_DIR / "histogram_ratio.png",
        save_stats=True,
        stats_path=stats_file,
    )
    plot_histogram(
        df,
        column="Seeds_per_Mbp",
        title="Distribution of Total Seeds per Mbp per Genome",
        xlabel="Total number of Seeds/Genome Size (Mbp) per Genome (log scale)",
        color="green",
        save_plot=True,
        plot_path=METABOLIC_POTENTIAL_PLOTS_DIR / "histogram_total_seeds_per_Mbp.png",
        save_stats=True,
        stats_path=stats_file,
        log_scale=True,
    )
    plot_histogram(
        df,
        column="Non_Seeds_per_Mbp",
        title="Distribution of Total Non-Seeds per Mbp per Genome",
        xlabel="Total number of non-Seeds/Genome Size (Mbp) per Genome (log scale)",
        color="lime",
        save_plot=True,
        plot_path=METABOLIC_POTENTIAL_PLOTS_DIR
        / "histogram_total_non_seeds_per_Mbp.png",
        save_stats=True,
        stats_path=stats_file,
        log_scale=True,
    )

    # Identify and plot outliers for seeds/non-seeds by phylum.
    outliers_ratio = get_outliers(
        df, column="Ratio", taxon="phylum", top_n=10, method="2std"
    )
    plot_outliers(
        df,
        outliers_ratio,
        save_plot=True,
        title=None,
        plot_path=METABOLIC_POTENTIAL_PLOTS_DIR / "outliers_phylum_ratio.png",
    )

    # Top/Bottom species by Ratio
    top_df, bottom_df = top_bottom_species(df, column="Ratio", top_n=10)
    top_df.to_csv(OUTPUT_DIR / "top_10_species_by_ratio.tsv", sep="\t", index=False)
    bottom_df.to_csv(
        OUTPUT_DIR / "bottom_10_species_by_ratio.tsv", sep="\t", index=False
    )
    plot_species_metrics(
        top_df,
        species_col="species",
        save_path=METABOLIC_POTENTIAL_PLOTS_DIR / "species_metrics_top_ratio.png",
        title="Metabolic Potential Metrics Across Top Species by Seed/Non-Seed Ratio",
    )
    plot_species_metrics(
        bottom_df,
        species_col="species",
        save_path=METABOLIC_POTENTIAL_PLOTS_DIR / "species_metrics_bottom_ratio.png",
        title="Metabolic Potential Metrics Across Bottom Species by Seed/Non-Seed Ratio",
    )

    # Plot violin plot for seeds/non-seeds ratio by phylum
    plot_violin_ratio_by_taxon(
        df,
        column="Ratio",
        top_n=10,
        y_limit=(0, 0.3),
        save_plot=True,
        plot_path=METABOLIC_POTENTIAL_PLOTS_DIR / "violin_ratio_phyla.png",
    )

    # Kruskal-Wallis test for Ratio by phylum (run once & log once)
    print("\nKruskal-Wallis test for Ratio per phylum:")
    H, p = kruskal_test_by_taxon(df, "Ratio", "phylum")
    with open(stats_file, "a", encoding="utf-8") as f:
        f.write(
            f"\n## Kruskal-Wallis (Ratio across phyla, top10)\nH={H:.3f}, p={p:.2e}\n"
        )

    # Dunn's post-hoc for Ratio by phylum
    dunn_out = OUTPUT_DIR / "dunn_posthoc_ratio_phylum.tsv"
    dunn_df = dunn_posthoc_test(
        df,
        taxon="phylum",
        column="Ratio",
        top_n=10,
        p_adjust="fdr_bh",
        export_tsv=dunn_out,
    )
    print("\nDunn's post-hoc test results (head):")
    print(dunn_df.head())

    # Dunn summary (write to TXT if there are results)
    if dunn_df.empty:
        with open(stats_file, "a", encoding="utf-8") as f:
            f.write("\n## Dunn post-hoc (Ratio across phyla, top10)\n")
            f.write("Not enough groups with data.\n")
    else:
        alpha = 0.05
        # Count only the upper triangle to avoid duplicates
        sig_mask = dunn_df.values < alpha
        n_sig = int(np.triu(sig_mask, k=1).sum())
        with open(stats_file, "a", encoding="utf-8") as f:
            f.write(f"\n## Dunn post-hoc (FDR-BH, α={alpha}, top10)\n")
            f.write(f"Significant pairs (unique): {n_sig}\n")
            f.write(f"Full p-value matrix saved to: {dunn_out}\n")

    # Heatmap of Dunn's test results
    plot_dunn_heatmap(
        dunn_df, save_path=METABOLIC_POTENTIAL_PLOTS_DIR / "dunn_posthoc_heatmap.png"
    )
    
    # Scatter: seeds vs non-seeds, colored by genome size
    plot_scatter_continuous_colormap(
        df,
        "Total_Seeds",
        "Total_non_Seeds",
        color_col="genome_size_Mbp",
        save_plot=True,
        plot_path=METABOLIC_POTENTIAL_PLOTS_DIR / "scatter_seed_vs_non_seed.png",
    )

    # Correlation stats (consistent with plot)
    pearson_r, p1, spearman_r, p2 = correlation_stats(
        df, x="Total_Seeds", y="Total_non_Seeds", log_x=False, log_y=False
        )

    with open(stats_file, "a", encoding="utf-8") as f:
        f.write("\n## Correlation (Seeds vs Non-Seeds)\n")
        f.write(f"Pearson r={pearson_r:.3f}, p={p1:.2e}\n")
        f.write(f"Spearman r={spearman_r:.3f}, p={p2:.2e}\n")

    # Scatter: Ratio vs genome size
    plot_scatter(
        df,
        title=None,
        y="Ratio",
        yscale="linear",
        save_plot=True,
        plot_path=METABOLIC_POTENTIAL_PLOTS_DIR / "scatter_ratio_vs_genome_size.png",
    )

    # Keep the top 10 phyla for stats as in plot
    top_taxa = df["phylum"].value_counts().nlargest(10).index
    df_sub = df[df["phylum"].isin(top_taxa)]

    # Correlation stats
    pearson_r, p1, spearman_r, p2 = correlation_stats(
        df_sub, x="genome_size_Mbp", y="Ratio", log_x=False, log_y=False
    )
    with open(stats_file, "a", encoding="utf-8") as f:
        f.write("\n## Correlation (Ratio vs genome_size_Mbp)\n")
        f.write(f"Pearson r={pearson_r:.3f}, p={p1:.2e}\n")
        f.write(f"Spearman r={spearman_r:.3f}, p={p2:.2e}\n")


if __name__ == "__main__":
    main()
