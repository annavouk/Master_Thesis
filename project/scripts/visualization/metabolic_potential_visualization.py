"""
Metabolic Potential Quantification (Approach 1) - Visualization

Generates histograms, boxplots, scatterplots, barplots, and violin plots to visualize
the metabolic potential (seeds/non-seeds/ratio) per genome and taxonomic group.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression   # <-- Εδώ το χρειαζόσουν!
from config import METABOLIC_POTENTIAL_1
from utils import load_data, split_and_clean_taxonomy
from scipy.stats import kruskal

def annotate_hist(ax, vals, label):
    median = vals.median()
    mean = vals.mean()
    std = vals.std()
    minv = vals.min()
    maxv = vals.max()
    ax.axvline(median, color='red', linestyle='--', label=f"Median = {median:.2f}")
    ax.axvline(mean, color='purple', linestyle=':', label=f"Mean = {mean:.2f}")
    ax.axvspan(mean-std, mean+std, color='purple', alpha=0.08, label=f"±1 STD = {std:.2f}")
    ax.legend()
    ax.annotate(f"Min: {minv:.2f}", xy=(minv, 0), xytext=(minv, 2), color='black', fontsize=8, rotation=90)
    ax.annotate(f"Max: {maxv:.2f}", xy=(maxv, 0), xytext=(maxv, 2), color='black', fontsize=8, rotation=90)

def annotate_box(ax, vals):
    median = vals.median()
    mean = vals.mean()
    std = vals.std()
    ax.text(0.98, 0.95,
        f"Mean: {mean:.2f}   Median: {median:.2f}   STD: {std:.2f}",
        color='black', fontsize=11, ha='right', va='top', transform=ax.transAxes)

def plot_distribution(df):
    """
    Plot frequency distributions and boxplot of seeds, non-seeds, and ratio per genome,
    with annotation for median, mean, std.
    """
    seeds = df['Total_Seeds'].dropna()
    nonseeds = df['Total_non_Seeds'].dropna()
    ratios = df['Ratio'].dropna()

    # Seeds Histogram + Annotation
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(seeds, bins=30, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Total number of Seeds per Genome')
    ax.set_ylabel('Number of Genomes')
    ax.set_title(f'Distribution of Total number of Seeds per Genome (N={len(seeds):,})')
    annotate_hist(ax, seeds, 'Seeds')
    plt.tight_layout()
    plt.show()

    # Non-seeds Histogram + Annotation
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(nonseeds, bins=30, edgecolor='black', alpha=0.7, color='skyblue')
    ax.set_xlabel('Total number of non-Seeds per Genome')
    ax.set_ylabel('Number of Genomes')
    ax.set_title(f'Distribution of Total number of non-Seeds per Genome (N={len(nonseeds):,})')
    annotate_hist(ax, nonseeds, 'Non-Seeds')
    plt.tight_layout()
    plt.show()

    # Ratio Histogram + Annotation
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(ratios, bins=30, edgecolor='black', alpha=0.7, color='orange')
    ax.set_xlabel('Seeds/Non-Seeds Ratio per Genome')
    ax.set_ylabel('Number of Genomes')
    ax.set_title(f'Distribution of Seeds/Non-Seeds Ratio per Genome (N={len(ratios):,})')
    annotate_hist(ax, ratios, 'Ratio')
    plt.tight_layout()
    plt.show()

    # Boxplot for Ratio + Annotation
    fig, ax = plt.subplots(figsize=(7, 4))
    sns.boxplot(x=ratios, color="lightgrey", ax=ax)
    ax.set_xlabel("Seeds/Non-Seeds Ratio")
    ax.set_title(f"Seeds/Non-Seeds Ratio per Genome (N={len(ratios):,})")
    annotate_box(ax, ratios)
    plt.tight_layout()
    plt.show()

    # Print summary statistics
    print("Total_Seeds summary:")
    print(seeds.describe())
    print("Total_non_Seeds summary:")
    print(nonseeds.describe())
    print("Ratio summary:")
    print(ratios.describe())

def analyze_ratio_std(df):
    mean = df['Ratio'].mean()
    std = df['Ratio'].std()
    lower = mean - std
    upper = mean + std

    below = df[(df['Ratio'] < lower)]
    above = df[(df['Ratio'] > upper)]

    print(f"Genomes BELOW mean-std: {len(below)} / {len(df)} ({100*len(below)/len(df):.2f}%)")
    print(f"Genomes ABOVE mean+std: {len(above)} / {len(df)} ({100*len(above)/len(df):.2f}%)")

    print("\nTop 10 phyla BELOW mean-std:")
    print(below['phylum'].value_counts().head(10))

    print("\nTop 10 phyla ABOVE mean+std:")
    print(above['phylum'].value_counts().head(10))

    return below, above 

def plot_extreme_phyla(below, above, top_n=10):
    below_counts = below['phylum'].value_counts().head(top_n)
    above_counts = above['phylum'].value_counts().head(top_n)

    fig, axs = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

    bars = axs[0].bar(below_counts.index, below_counts.values, color='blue', edgecolor='black')
    axs[0].set_title(f"Phyla BELOW mean-std\n(N={len(below)})")
    axs[0].set_xlabel("Phylum")
    axs[0].set_ylabel("Number of Genomes")
    axs[0].tick_params(axis='x', rotation=45)
    for bar in bars:
        height = bar.get_height()
        axs[0].annotate(f'{int(height)}', xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=10)

    bars = axs[1].bar(above_counts.index, above_counts.values, color='red', edgecolor='black')
    axs[1].set_title(f"Phyla ABOVE mean+std\n(N={len(above)})")
    axs[1].set_xlabel("Phylum")
    axs[1].tick_params(axis='x', rotation=45)
    for bar in bars:
        height = bar.get_height()
        axs[1].annotate(f'{int(height)}', xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    plt.show()

def plot_scatter_ratio_vs_genome(df):
    """
    Scatter plot of Ratio vs genome size by phylum (top 10 only for color clarity).
    """
    top_phyla = df['phylum'].value_counts().nlargest(10).index
    sub = df[df['phylum'].isin(top_phyla)]
    n = len(sub)
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=sub, x='genome_length', y='Ratio', hue='phylum', alpha=0.7)
    plt.title(f'Variation in Seeds/Non-Seeds Ratio with Genome Size among top 10 Phyla\n(N = {n} genomes)')
    plt.xlabel('Genome Size (bp)')
    plt.ylabel('Seed/non-Seed Ratio')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Phylum')
    plt.tight_layout()
    plt.show()


def plot_scatter_nonseeds_vs_genome(df):
    """
    Scatter plot of Total_non_Seeds vs genome size by phylum (top 10 only for color clarity).
    """
    top_phyla = df['phylum'].value_counts().nlargest(10).index
    sub = df[df['phylum'].isin(top_phyla)]
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=sub, x='genome_length', y='Total_non_Seeds', hue='phylum', alpha=0.7)
    plt.title('Total number of non-Seeds vs Genome Size by Phylum')
    plt.xlabel('Genome Size (bp)')
    plt.ylabel('Total number of non-Seeds')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title='Phylum')
    plt.tight_layout()
    plt.show()

def plot_seeds_vs_nonseeds_with_regression(df):
    """
    Scatter plot Seeds vs Non-Seeds, χρώμα genome size, regression line, correlation annotations.
    """
    plt.figure(figsize=(7, 5))
    sc = plt.scatter(
        df['Total_Seeds'], df['Total_non_Seeds'],
        c=df['genome_length'], cmap='viridis', s=10, alpha=0.5
    )

    # Regression line (OLS)
    x = df['Total_Seeds'].values.reshape(-1, 1)
    y = df['Total_non_Seeds'].values
    model = LinearRegression().fit(x, y)
    x_line = np.linspace(x.min(), x.max(), 200).reshape(-1, 1)
    y_line = model.predict(x_line)
    plt.plot(x_line, y_line, color='red', lw=2, label='Regression line')

    # Pearson/Spearman r
    pearson_corr, p_pearson = pearsonr(df['Total_Seeds'], df['Total_non_Seeds'])
    spearman_corr, p_spearman = spearmanr(df['Total_Seeds'], df['Total_non_Seeds'])

    textstr = (
        f"Pearson r = {pearson_corr:.2f} (p = {p_pearson:.2e})\n"
        f"Spearman r = {spearman_corr:.2f} (p = {p_spearman:.2e})"
    )
    plt.gca().text(
        0.05, 0.95, textstr,
        transform=plt.gca().transAxes,
        fontsize=11, verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="gray", alpha=0.8)
    )

    plt.xlabel('Number of Seeds per Genome')
    plt.ylabel('Number of Non-Seeds per Genome')
    plt.title('Seeds vs Non-Seeds per Genome (colored by Genome Size)')
    cbar = plt.colorbar(sc)
    cbar.set_label('Genome Size (bp)')
    plt.legend()
    plt.tight_layout()
    plt.show()

def regression_nonseeds_on_seeds_and_genomesize(df):
    X = df[['Total_Seeds', 'genome_length']]
    X = sm.add_constant(X)
    y = df['Total_non_Seeds']

    model = sm.OLS(y, X).fit()
    print(model.summary())
    return model

def plot_top_genera_barplots(df):
    """
    Generate barplots for top 10 genera with highest/lowest number of seeds and non-seeds.
    """
    max_seeds = df.groupby('genus')['Total_Seeds'].max().sort_values(ascending=False).head(10)
    min_seeds = df.groupby('genus')['Total_Seeds'].min().sort_values().head(10)
    max_non_seeds = df.groupby('genus')['Total_non_Seeds'].max().sort_values(ascending=False).head(10)
    min_non_seeds = df.groupby('genus')['Total_non_Seeds'].min().sort_values().head(10)

    plt.figure(figsize=(10, 5))
    max_seeds.plot(kind='bar', color='green', edgecolor='black')
    plt.title("Top 10 Genera with Highest Number of Seeds in a Single Genome")
    plt.ylabel("Number of Seeds")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 5))
    min_seeds.plot(kind='bar', color='limegreen', edgecolor='black')
    plt.title("Top 10 Genera with Lowest Number of Seeds in a Single Genome")
    plt.ylabel("Number of Seeds")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 5))
    max_non_seeds.plot(kind='bar', color='orange', edgecolor='black')
    plt.title("Top 10 Genera with Highest Number of non-Seeds in a Single Genome")
    plt.ylabel("Number of non-Seeds")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 5))
    min_non_seeds.plot(kind='bar', color='salmon', edgecolor='black')
    plt.title("Top 10 Genera with Lowest Number of non-Seeds in a Single Genome")
    plt.ylabel("Number of non-Seeds")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()


def plot_outlier_taxa_bars(df, top_n=10):
    """
    Plot the top 10 genera/classes in the top/bottom 5% for seeds and non-seeds.
    """
    def plot_top_taxa(counts, title, color):
        plt.figure(figsize=(10, 6))
        counts.head(10).plot(kind='bar', color=color, edgecolor='black')
        plt.title(title)
        plt.ylabel('Count')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.show()

    lower_seeds = df['Total_Seeds'].quantile(0.05)
    upper_seeds = df['Total_Seeds'].quantile(0.95)
    lower_nonseeds = df['Total_non_Seeds'].quantile(0.05)
    upper_nonseeds = df['Total_non_Seeds'].quantile(0.95)

    # Seeds
    low_seeds = df[df['Total_Seeds'] <= lower_seeds]
    high_seeds = df[df['Total_Seeds'] >= upper_seeds]

    plot_top_taxa(low_seeds['genus'].value_counts().nlargest(top_n), 'Top 10 Genera in Low Seeds Outliers', 'lightcoral')
    plot_top_taxa(low_seeds['class'].value_counts().nlargest(top_n), 'Top 10 Classes in Low Seeds Outliers', 'skyblue')
    plot_top_taxa(high_seeds['genus'].value_counts().nlargest(top_n), 'Top 10 Genera in High Seeds Outliers', 'lightcoral')
    plot_top_taxa(high_seeds['class'].value_counts().nlargest(top_n), 'Top 10 Classes in High Seeds Outliers', 'skyblue')

    # Non-Seeds
    low_nonseeds = df[df['Total_non_Seeds'] <= lower_nonseeds]
    high_nonseeds = df[df['Total_non_Seeds'] >= upper_nonseeds]

    plot_top_taxa(low_nonseeds['genus'].value_counts().nlargest(top_n), 'Top 10 Genera in Low Non-Seeds Outliers', 'lightcoral')
    plot_top_taxa(low_nonseeds['class'].value_counts().nlargest(top_n), 'Top 10 Classes in Low Non-Seeds Outliers', 'skyblue')
    plot_top_taxa(high_nonseeds['genus'].value_counts().nlargest(top_n), 'Top 10 Genera in High Non-Seeds Outliers', 'lightcoral')
    plot_top_taxa(high_nonseeds['class'].value_counts().nlargest(top_n), 'Top 10 Classes in High Non-Seeds Outliers', 'skyblue')


def plot_violin_ratio_by_phylum(df, top_n=10):
    """
    Generate violin plots for total number of seeds to non-seeds per phylum.
    """
    top_taxa = df['phylum'].value_counts().nlargest(top_n).index
    df_sub = df[df['phylum'].isin(top_taxa)]

    n_per_phylum = df_sub['phylum'].value_counts().reindex(top_taxa)
    avg_size_per_phylum = df_sub.groupby('phylum')['genome_length'].mean().reindex(top_taxa).round(0).astype(int)

    medians = df_sub.groupby('phylum')['Ratio'].median().reindex(top_taxa)

    sns.set(style="whitegrid")

    plt.figure(figsize=(14, 6))
    ax = sns.violinplot(data=df_sub, x='phylum', y='Ratio', inner='box', palette='pastel', legend=False, order=top_taxa)
    plt.title('Distribution of Seed/non-Seed Ratio per Phylum')
    plt.xticks(rotation=45, ha='right')

    xticklabels = [
        f"{taxon}\n(N={n_per_phylum[taxon]:,}, avg size={avg_size_per_phylum[taxon]:,} bp)"
        for taxon in top_taxa
    ]
    ax.set_xticklabels(xticklabels)

    for i, phylum in enumerate(top_taxa):
        ax.scatter(i, medians[phylum], color='black', marker='o', s=50, zorder=10, label='Median' if i==0 else "")

    plt.tight_layout()
    plt.show()

def kruskal_test_by_group(df, value_col, group_col, top_n=10):
    top_groups = df[group_col].value_counts().nlargest(top_n).index
    sub = df[df[group_col].isin(top_groups)]
    data = [sub[sub[group_col] == group][value_col].dropna() for group in top_groups]
    H, p = kruskal(*data)
    print(f"Kruskal–Wallis for '{value_col}' by '{group_col}': H={H:.2f}, p={p:.3g}")
    return H, p


def main():
    # Load and split taxonomy
    df = load_data(METABOLIC_POTENTIAL_1, filetype="csv")
    taxonomy_df = split_and_clean_taxonomy(df, 'gtdb_taxonomy')
    df = pd.concat([df, taxonomy_df], axis=1)

    #plot_distribution(df)

    below, above = analyze_ratio_std(df)
    #plot_extreme_phyla(below, above, top_n=10)

    phyla_above = above['phylum'].value_counts().head(10).index
    phyla_below = below['phylum'].value_counts().head(10).index
    outlier_phyla = set(phyla_above) | set(phyla_below)

    total_counts = df['phylum'].value_counts()

    #print("Phylum\tTotal N\tAbove mean+std\tBelow mean-std")
    #for phylum in outlier_phyla:
    #    total = total_counts.get(phylum, 0)
    #    above_n = above[above['phylum'] == phylum].shape[0]
    #    below_n = below[below['phylum'] == phylum].shape[0]
    #    print(f"{phylum}\t{total}\t{above_n}\t{below_n}")

    #plot_scatter_ratio_vs_genome(df)

    plot_seeds_vs_nonseeds_with_regression(df)
    regression_nonseeds_on_seeds_and_genomesize(df)

    #Will be deleted 
    #plot_top_genera_barplots(df)
    #plot_outlier_taxa_bars(df)

    plot_violin_ratio_by_phylum(df, top_n=10)

    # Kruskal-Wallis test for Ratio by phylum
    print("\nKruskal–Wallis test for Ratio per phylum:")
    kruskal_test_by_group(df, 'Ratio', 'phylum')
    

if __name__ == "__main__":
    main()
