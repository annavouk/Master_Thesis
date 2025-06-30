"""

"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from config import SEEDS_PICKLE, NON_SEEDS_PICKLE, AMINOACIDS, COMPACT_METADATA, OUTPUT_DIR
from utils import load_data, split_and_clean_taxonomy


def main():
    # Load data
    seeds = load_data(SEEDS_PICKLE, filetype='pickle')
    non_seeds = load_data(NON_SEEDS_PICKLE, filetype='pickle')
    aa_df = load_data(AMINOACIDS, filetype='tsv')
    metadata_df = load_data(COMPACT_METADATA, filetype='csv')
    
    aa_seed_ids = aa_df["SEED_ID"].tolist()
   
    aa_in_seed_matrix = [cpd for cpd in aa_seed_ids if cpd in seeds.columns]
    seed_aa_matrix = seeds[aa_in_seed_matrix]

    aa_in_non_seed_matrix = [cpd for cpd in aa_seed_ids if cpd in non_seeds.columns]
    nonseed_aa_matrix = non_seeds[aa_in_non_seed_matrix]

    # counts per amino acid for seeds and non-seeds
    seed_counts = seed_aa_matrix.sum(axis=0)
    nonseed_counts = nonseed_aa_matrix.sum(axis=0)

    # Combine into DataFrame for grouped barplot
    counts_df = pd.DataFrame({
        'Auxotrophy (Seed)': seed_counts,
        'Prototrophy (Non-Seed)': nonseed_counts
    })


    labels = aa_df.set_index('SEED_ID')['compound_name'].reindex(counts_df.index)
    counts_df.index = labels

    # Bar Plot: Auxotrophy/Prototrophy Frequency per Amino Acid
    counts_df.plot(kind='bar', figsize=(12, 6))
    plt.title('Auxotrophy/Prototrophy Frequency per Amino Acid')
    plt.ylabel('Number of Genomes')
    plt.xlabel('Amino Acid (SEED_ID)')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.show()

    # Count auxotrophy/prototrophy per genome
    auxotrophy_per_genome = seed_aa_matrix.sum(axis=1)
    prototrophy_per_genome = nonseed_aa_matrix.sum(axis=1)

    summary_df = pd.DataFrame({
        'Auxotrophies': auxotrophy_per_genome,
        'Prototrophies': prototrophy_per_genome
    })

    #print(summary_df.describe())

    # Histogram: Auxotrophies and Prototrophies per Genome
    summary_df.plot(kind='hist', bins=range(0, len(aa_seed_ids)+2), alpha=0.7)
    plt.xlabel('Number of Amino Acids')
    plt.ylabel('Number of Genomes')
    plt.title('Distribution of Auxotrophies and Prototrophies per Genome')
    plt.tight_layout()
    plt.show()

    # Genomes with most seed amino acids (auxotrophies)
    top_auxotrophy = auxotrophy_per_genome.sort_values(ascending=False).head(10)
    #print("Top genomes with most auxotrophies (seeds):")
    #print(top_auxotrophy)

    # Genomes with most non-seed amino acids (prototrophies)
    top_prototrophy = prototrophy_per_genome.sort_values(ascending=False).head(10)
    #print("\nTop genomes with most prototrophies (non-seeds):")
    #print(top_prototrophy)

    # Barplot: Top genomes (auxotrophies)
    top_auxotrophy.plot(kind='bar')
    plt.ylabel('Number of Auxotrophies (seed amino acids)')
    plt.xlabel('Genome')
    plt.title('Top 10 Genomes with Most Auxotrophies')
    plt.tight_layout()
    plt.show()

    # Merge
    seed_aa_matrix.index = seed_aa_matrix.index.astype(str)
    metadata_df['patric_id'] = metadata_df['patric_id'].astype(str)
    merged = seed_aa_matrix.merge(metadata_df[['patric_id', 'gtdb_taxonomy']], 
                             left_index=True, right_on='patric_id', how='left')
    
    # Split and clean taxonomy
    taxonomy_df = split_and_clean_taxonomy(merged, 'gtdb_taxonomy')
    merged = pd.concat([merged, taxonomy_df], axis=1)

    # Count amino acids as seeds per phylum
    phylum_aa_summary = merged.groupby('phylum')[aa_in_seed_matrix].sum()
    #print(phylum_aa_summary.head())

    # Calculate percentage of auxotrophies per phylum
    phylum_counts = merged['phylum'].value_counts()
    phylum_aa_percentage = phylum_aa_summary.div(phylum_counts, axis=0) * 100
    #print(phylum_aa_percentage.head())

    # Stacked bar Plot: Auxotrophies per Phylum
    ax = phylum_aa_summary.T.plot(kind='bar', stacked=True, figsize=(15,7))
    plt.ylabel('Number of Genomes')
    plt.title('Auxotrophy (Seed) Frequency per Amino Acid per Phylum')
    plt.tight_layout()
    plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=4, fontsize=8)
    plt.show()


if __name__ == "__main__":
    main()