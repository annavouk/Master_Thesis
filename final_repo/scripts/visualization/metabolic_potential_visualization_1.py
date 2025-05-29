"""
Visualize the metabolic potential of each genome.
Generate plots...
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def split_and_clean_taxonomy(df, taxonomy_col):
    """
    Split gtdb_taxonomy string to taxonomic levels.

    Args:
        df (pd.Dataframe): gtdb_taxonomy column with each level started with e.g. d__
        for domain and seperated by ";".

    Returns:
        df (pd.DataFrame): A new DataFrame where each taxonomy rank occupies its own column,
        with prefixes removed.
    """
    taxonomy_series = df.loc[:, taxonomy_col].astype(str)
    taxonomy_split = taxonomy_series.str.split(';', expand=True)

    levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxonomy_split.columns = levels[:taxonomy_split.shape[1]]

    taxonomy_cleaned = taxonomy_split.apply(lambda col: col.str.replace(r'^[a-z]__', '', regex=True))
    return taxonomy_cleaned


if __name__ == "__main__":

    # Load the data
    df = pd.read_csv("../analysis/metabolic_potential_summary.csv", low_memory=False)


    """
    Generate histograms and boxplot to visualize the distribution of Seeds,non-Seeds and Ratio per genome.
    """
    # Group by number of seeds and non-seeds and count genomes
    seed_counts = df['Total_Seeds'].value_counts().sort_index()
    non_seed_counts = df['Total_non_Seeds'].value_counts().sort_index()

    # Plot side-by-side
    fig, axs = plt.subplots(1, 2, figsize=(14, 12), sharey=True)

    # Seeds histogram plot
    axs[0].bar(seed_counts.index, seed_counts.values, color='green', edgecolor='black')
    axs[0].set_xlabel('Number of Seeds')
    axs[0].set_ylabel('Number of Genomes')
    axs[0].set_title('Distribution of Seeds per Genome')

    # Non-Seeds histogram plot
    axs[1].bar(non_seed_counts.index, non_seed_counts.values, color='skyblue', edgecolor='black')
    axs[1].set_xlabel('Number of non-Seeds')
    axs[1].set_title('Distribution of non-Seeds per Genome')

    plt.tight_layout()
    plt.show()

    # Boxplot Ratio per Genome
    plt.figure(figsize=(6, 4))
    sns.boxplot(x=df["Ratio"])
    plt.xlabel("Seeds/non-Seeds Ratio")
    plt.title("Seeds/non-Seeds Ratio per Genome")
    plt.show()

    # Print summary statistics
    #print(df["Total_Seeds"].describe())
    #print(df["Total_non_Seeds"].describe())
    #print(df["Ratio"].describe())


    """
    Use the split and clean taxonomy function to rearange the df and merge taxonomy back to df.
    """
    taxonomy_df = split_and_clean_taxonomy(df, 'gtdb_taxonomy')

    df = pd.concat([df, taxonomy_df], axis=1)


    """
    Generate scatter plots to visualize the correlation of total number of seeds/non-seeds with genome size
    and taxonomic group.
    """
    # Filter only top 10 phyla for clarity
    top_phyla = df['phylum'].value_counts().nlargest(10).index
    df = df[df['phylum'].isin(top_phyla)]

    # Scatter plot - Seeds
    plt.figure(figsize=(10, 6))
    #sns.scatterplot(data=df, x='genome_length', y='Total_Seeds', hue='phylum', alpha=0.7)
    sns.scatterplot(data=df, x='genome_length', y='Total_non_Seeds', hue='phylum', alpha=0.7)
    #plt.title('Total number of Seeds vs Genome Size by Phylum')
    plt.title('Total number of non-Seeds vs Genome Size by Phylum')
    plt.xlabel('Genome Size (bp)')
    #plt.ylabel('Total number of Seeds')
    plt.ylabel('Total number of non-Seeds')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()

    """
    Generate bar plot to visualize the top 10 most abundand genera collecting the highest/lowest number of Seeds/non-Seeds.
    """
    # Top 10 genera by max seeds
    max_seeds = df.groupby('genus')['Total_Seeds'].max().sort_values(ascending=False).head(10)

    # Top 10 genera by min seeds
    min_seeds = df.groupby('genus')['Total_Seeds'].min().sort_values().head(10)

    # Top 10 genera by max non-seeds
    max_non_seeds = df.groupby('genus')['Total_non_Seeds'].max().sort_values(ascending=False).head(10)

    # Top 10 genera by min non-seeds
    min_non_seeds = df.groupby('genus')['Total_non_Seeds'].min().sort_values().head(10)

    # Plot: Top 10 genera with highest number of seeds
    plt.figure(figsize=(10, 5))
    max_seeds.plot(kind='bar', color='green', edgecolor='black')
    plt.title("Top 10 Genera with Highest Number of Seeds in a Single Genome")
    plt.ylabel("Number of Seeds")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    # Plot: Top 10 genera with lowest number of seeds
    plt.figure(figsize=(10, 5))
    min_seeds.plot(kind='bar', color='limegreen', edgecolor='black')
    plt.title("Top 10 Genera with Lowest Number of Seeds in a Single Genome")
    plt.ylabel("Number of Seeds")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    # Plot: Top 10 genera with highest number of non-seeds
    plt.figure(figsize=(10, 5))
    max_non_seeds.plot(kind='bar', color='orange', edgecolor='black')
    plt.title("Top 10 Genera with Highest Number of non-Seeds in a Single Genome")
    plt.ylabel("Number of non-Seeds")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    # Plot: Top 10 genera with lowest number of non-seeds
    plt.figure(figsize=(10, 5))
    min_non_seeds.plot(kind='bar', color='salmon', edgecolor='black')
    plt.title("Top 10 Genera with Lowest Number of non-Seeds in a Single Genome")
    plt.ylabel("Number of non-Seeds")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    """
    Based on the histogram of distribution of total number of seeds per genome,
    explore the 5% of the sample collecting the highest/lowest number of seeds/non-seeds.
    """
    # Calculate thresholds
    lower_threshold = df['Total_Seeds'].quantile(0.05)
    upper_threshold = df['Total_Seeds'].quantile(0.95)

    # Calculate thresholds
    lower_threshold_non_seeds = df['Total_non_Seeds'].quantile(0.05)
    upper_threshold_non_seeds = df['Total_non_Seeds'].quantile(0.95)

    # Subset for bottom 5% (low outliers) and top 5% (high outliers) based on seeds count
    low_outliers = df[df['Total_Seeds'] <= lower_threshold]
    high_outliers = df[df['Total_Seeds'] >= upper_threshold]

    # Subset for bottom 5% (low outliers) and top 5% (high outliers) based on non-seeds count
    low_outliers_non_seeds = df[df['Total_non_Seeds'] <= lower_threshold_non_seeds]
    high_outliers_non_seeds = df[df['Total_non_Seeds'] >= upper_threshold_non_seeds]

    # For low outliers (e.g., bottom 5% by total non-seeds)
    #low_genera_counts = low_outliers_non_seeds['genus'].value_counts()
    #low_classes_counts = low_outliers_non_seeds['class'].value_counts()

    #print(f"Low outliers - Unique genera: {low_outliers_non_seeds['genus'].nunique()}")
    #print(low_genera_counts.head(10))

    #print(f"Low outliers - Unique classes: {low_outliers_non_seeds['class'].nunique()}")
    #print(low_classes_counts.head(10))


    # For high outliers (e.g., top 5% by total non-seeds)
    #high_genera_counts = high_outliers_non_seeds['genus'].value_counts()
    #high_classes_counts = high_outliers_non_seeds['class'].value_counts()

    #print(f"High outliers - Unique genera: {high_outliers_non_seeds['genus'].nunique()}")
    #print(high_genera_counts.head(10))

    #print(f"High outliers - Unique classes: {high_outliers_non_seeds['class'].nunique()}")
    #print(high_classes_counts.head(10))

    # Plot
    def plot_top_taxa(counts, title, color):
        plt.figure(figsize=(10, 6))
        counts.head(10).plot(kind='bar', color=color, edgecolor='black')
        plt.title(title)
        plt.ylabel('Count')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.show()


    # --- Seeds ---
    # Low outliers seeds
    low_seeds_genera = low_outliers['genus'].value_counts()
    low_seeds_classes = low_outliers['class'].value_counts()

    plot_top_taxa(low_seeds_genera, 'Top 10 Genera in Low Seeds Outliers', 'lightcoral')
    plot_top_taxa(low_seeds_classes, 'Top 10 Classes in Low Seeds Outliers', 'skyblue')

    # High outliers seeds
    high_seeds_genera = high_outliers['genus'].value_counts()
    high_seeds_classes = high_outliers['class'].value_counts()

    plot_top_taxa(high_seeds_genera, 'Top 10 Genera in High Seeds Outliers', 'lightcoral')
    plot_top_taxa(high_seeds_classes, 'Top 10 Classes in High Seeds Outliers', 'skyblue')


    # --- Non-Seeds ---
    # Low outliers non-seeds
    low_non_seeds_genera = low_outliers_non_seeds['genus'].value_counts()
    low_non_seeds_classes = low_outliers_non_seeds['class'].value_counts()

    plot_top_taxa(low_non_seeds_genera, 'Top 10 Genera in Low Non-Seeds Outliers', 'lightcoral')
    plot_top_taxa(low_non_seeds_classes, 'Top 10 Classes in Low Non-Seeds Outliers', 'skyblue')

    # High outliers non-seeds
    high_non_seeds_genera = high_outliers_non_seeds['genus'].value_counts()
    high_non_seeds_classes = high_outliers_non_seeds['class'].value_counts()

    plot_top_taxa(high_non_seeds_genera, 'Top 10 Genera in High Non-Seeds Outliers', 'lightcoral')
    plot_top_taxa(high_non_seeds_classes, 'Top 10 Classes in High Non-Seeds Outliers', 'skyblue')

    """
    Generate violin plots to visualize number of seeds/non-seeds per phyla.
    """
    # Set plot style
    sns.set(style="whitegrid")

    # Violin plot for Total_Seeds per phylum
    plt.figure(figsize=(14, 6))
    sns.violinplot(data=df, x='phylum', y='Total_Seeds', inner='box', palette='pastel')
    plt.title('Distribution of Total Seeds per Phylum')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    # Violin plot for Total_non_Seeds per phylum
    plt.figure(figsize=(14, 6))
    sns.violinplot(data=df, x='phylum', y='Total_non_Seeds', inner='box', palette='muted')
    plt.title('Distribution of Total Non-Seeds per Phylum')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()
