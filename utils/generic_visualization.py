import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


# ------------------------
# Histogram
# ------------------------
def annotate_hist(ax, vals, label, method="2std", log_scale=False):
    """
    Annotate histogram with statistics and optional outlier shading.
    If log_scale=True, compute stats in log10 space.
    """
    n = len(vals)

    data = np.log10(vals[vals > 0]) if log_scale else vals
    median = data.median()
    mean = data.mean()
    std = data.std()
    minv = data.min()
    maxv = data.max()

    ax.annotate(
        f"N = {n:,} genomes",
        xy=(0.99, 0.9), xycoords="axes fraction",
        fontsize=8, color="grey", ha="right", va="center"
    )

    # Median & Mean lines
    ax.axvline(10**median if log_scale else median,
               color="red", linestyle="--", label=f"Median = {median:.2f}")
    ax.axvline(10**mean if log_scale else mean,
               color="purple", linestyle=":", label=f"Mean = {mean:.2f}")

    # Shading depending on method
    if method == "2std":
        lower, upper = mean - 2*std, mean + 2*std
        if log_scale:
            lower, upper = 10**lower, 10**upper
        ax.axvspan(lower, upper, color="purple", alpha=0.08,
                   label=f"±2 STD (σ={std:.2f})")
    elif method == "iqr":
        Q1, Q3 = data.quantile(0.25), data.quantile(0.75)
        if log_scale:
            Q1, Q3 = 10**Q1, 10**Q3
        ax.axvspan(Q1, Q3, color="green", alpha=0.08,
                   label=f"IQR [{Q1:.2f}, {Q3:.2f}]")

    # Min/Max
    ax.annotate(f"Min: {minv:.2f}", xy=(minv, 0), xytext=(minv, 15),
                color="black", fontsize=8, rotation=90)
    ax.annotate(f"Max: {maxv:.2f}", xy=(maxv, 0), xytext=(maxv, 15),
                color="black", fontsize=8, rotation=90)

    ax.legend()


def plot_histogram(
    df,
    column,
    bins=30,
    color="teal",
    title=None,
    xlabel=None,
    save_plot=False,
    plot_path=None,
    save_stats=False,
    stats_path=None,
    log_scale=False,
    xlim=None,
    method="2std",
):
    """Plot histogram of a specified column in a DataFrame."""
    vals = df[column].dropna()
    if log_scale:
        vals = vals[vals > 0]
        bins = np.logspace(np.log10(vals.min()), np.log10(vals.max()), bins)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(vals, bins=bins, edgecolor="black", alpha=0.7, color=color)
    if log_scale:
        ax.set_xscale("log")
        mean = vals.mean()
        std = vals.std()
        upper_limit = mean + 3 * std
        if xlim is None:
            xlim = (vals.min(), min(vals.max(), upper_limit))
        ax.set_xlim(xlim)
        if vals.max() > xlim[1]:
            ax.annotate(
                f"Max: {vals.max():.2f}",
                xy=(0.99, 0.95),
                xycoords="axes fraction",
                ha="right",
                va="top",
                fontsize=8,
                color="black",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7),
            )
    ax.set_xlabel(xlabel if xlabel else column)
    ax.set_ylabel("Number of Genomes")
    ax.set_title(
        title
        if title
        else f'Distribution of {column}{" (log scale)" if log_scale else ""} (N={len(vals):,} genomes)'
    )
    annotate_hist(ax, vals, column, method=method)
    plt.tight_layout()
    if save_plot and plot_path:
        plt.savefig(plot_path, dpi=300)
    plt.show()

    stats = vals.describe()
    print(f"{column} summary{' (log scale)' if log_scale else ''}:\n{stats}")
    if save_stats and stats_path:
        with open(stats_path, "a") as f:
            f.write(
                f"Summary statistics for {column}{' (log scale)' if log_scale else ''}:\n{stats}\n\n"
            )
    return stats


# ------------------------
# Boxplot
# ------------------------
def annotate_box(ax, vals, log_scale=False, method="2std"):
    """
    Annotate boxplot with mean/median/std or IQR shading.
    """
    data = np.log10(vals[vals > 0]) if log_scale else vals
    median = data.median()
    mean = data.mean()
    std = data.std()

    ax.text(0.95, 0.95,
            f"Mean: {mean:.2f}\nMedian: {median:.2f}\nSTD: {std:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="grey", alpha=0.6))

    if method == "2std":
        lower, upper = mean - 2*std, mean + 2*std
        if log_scale:
            lower, upper = 10**lower, 10**upper
        ax.axvspan(lower, upper, color="purple", alpha=0.08,
                   label=f"±2 STD (σ={std:.2f})", zorder=0)
    elif method == "iqr":
        Q1, Q3 = data.quantile(0.25), data.quantile(0.75)
        if log_scale:
            Q1, Q3 = 10**Q1, 10**Q3
        ax.axvspan(Q1, Q3, color="green", alpha=0.08,
                   label=f"IQR [{Q1:.2f}, {Q3:.2f}]", zorder=0)

    ax.legend(loc="upper left", fontsize=8)


def plot_boxplot(
    df,
    column,
    color="lightsteelblue",
    title=None,
    xlabel=None,
    ylabel=None,
    save_plot=False,
    plot_path=None,
    orient="h",
    log_scale=False,
    method="2std",
):
    """Plot boxplot of a specified column in a DataFrame."""
    vals = df[column].dropna()
    if log_scale:
        vals = np.log10(vals + 1e-9)
        xlabel = xlabel or f"Log10({column})"
    else:
        xlabel = xlabel or column

    fig, ax = plt.subplots(figsize=(7, 4))
    sns.boxplot(
        x=vals if orient == "h" else None,
        y=vals if orient == "v" else None,
        color=color,
        ax=ax,
    )
    ax.set_xlabel(xlabel if xlabel else (column if orient == "h" else ""))
    ax.set_ylabel(ylabel if ylabel else (column if orient == "v" else ""))
    ax.set_title(title if title else f"{column} (N={len(vals):,})")

    annotate_box(ax, vals, log_scale=log_scale, method=method)

    plt.tight_layout()
    if save_plot and plot_path:
        plt.savefig(plot_path, dpi=300)
    plt.show()
    return vals.describe()
