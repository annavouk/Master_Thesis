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
from config import METABOLIC_POTENTIAL_1
from utils import load_data, split_and_clean_taxonomy
from scipy.stats import kruskal


def plot_distribution(df):
    """
    Plot frequency distribution of Total_Seeds and Total_non_Seeds per genome,
    and a boxplot of the Ratio.

    Args:
        df (pd.DataFrame): Contains 'Total_Seeds', 'Total_non_Seeds', 'Ratio'.
    """
    seed_counts = df['Total_Seeds'].value_counts().sort_index()
    non_seed_counts = df['Total_non_Seeds'].value_counts().sort_index()

    fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    axs[0].bar(seed_counts.index, seed_counts.values, color='green', edgecolor='black')
    axs[0].set_xlabel('Number of Seeds')
    axs[0].set_ylabel('Number of Genomes')
    axs[0].set_title('Distribution of Seeds per Genome')

    axs[1].bar(non_seed_counts.index, non_seed_counts.values, color='skyblue', edgecolor='black')
    axs[1].set_xlabel('Number of non-Seeds')
    axs[1].set_title('Distribution of non-Seeds per Genome')

    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(6, 4))
    sns.boxplot(x=df["Ratio"], color="lightgrey")
    plt.xlabel("Seeds/non-Seeds Ratio")
    plt.title("Seeds/non-Seeds Ratio per Genome")
    plt.tight_layout()
    plt.show()

    # Print summary statistics
    print("Total_Seeds summary:")
    print(df["Total_Seeds"].describe())
    print("Total_non_Seeds summary:")
    print(df["Total_non_Seeds"].describe())
    print("Ratio summary:")
    print(df["Ratio"].describe())


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


def plot_violin_seed_nonseed_by_phylum(df, top_n=10):
    """
    Generate violin plots for total number of seeds and non-seeds per phylum.
    """
    top_taxa = df['phylum'].value_counts().nlargest(top_n).index
    df_sub = df[df['phylum'].isin(top_taxa)]

    sns.set(style="whitegrid")

    plt.figure(figsize=(14, 6))
    sns.violinplot(data=df_sub, x='phylum', y='Total_Seeds', inner='box', palette='pastel', hue='phylum', legend=False)
    plt.title('Distribution of Total Seeds per Phylum')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(14, 6))
    sns.violinplot(data=df_sub, x='phylum', y='Total_non_Seeds', inner='box', palette='muted', hue='phylum', legend=False)
    plt.title('Distribution of Total Non-Seeds per Phylum')
    plt.xticks(rotation=45, ha='right')
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

    plot_distribution(df)
    plot_scatter_nonseeds_vs_genome(df)
    plot_top_genera_barplots(df)
    plot_outlier_taxa_bars(df)
    plot_violin_seed_nonseed_by_phylum(df, top_n=10)

    # Kruskal-Wallis test for Total_Seeds and Total_non_Seeds by phylum
    print("\nKruskal–Wallis test for Total_Seeds per phylum:")
    kruskal_test_by_group(df, 'Total_Seeds', 'phylum')
    print("\nKruskal–Wallis test for Total_non_Seeds per phylum:")  
    kruskal_test_by_group(df, 'Total_non_Seeds', 'phylum')
    

if __name__ == "__main__":
    main()
