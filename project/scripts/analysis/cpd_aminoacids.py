"""
Analyze amino acid auxotrophy and prototrophy profiles across genomes.

Generate plots:
    - Barplot of auxotrophy/prototrophy frequency per amino acid
    - Histogram of auxotrophies/prototrophies per genome
    - Stacked barplots of auxotrophy patterns per taxonomic rank
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from config import SEEDS_PICKLE, NON_SEEDS_PICKLE, COMPACT_METADATA_WITH_BIOME, AMINOACIDS, OUTPUT_DIR
from utils import load_data, split_and_clean_taxonomy


# ------------------------
# Auxotrophy analysis
# ------------------------
def analyze_auxotrophy(seed_matrix, nonseed_matrix, aa_df, metadata_df):
    """Compute auxotrophy/prototrophy distributions and taxonomic aggregations."""
    aa_seed_ids = aa_df["SEED_ID"].tolist()

    # Filter matrices to amino acid columns
    aa_in_seed_matrix = [cpd for cpd in aa_seed_ids if cpd in seed_matrix.columns]
    aa_in_non_seed_matrix = [cpd for cpd in aa_seed_ids if cpd in nonseed_matrix.columns]
    seed_aa_matrix = seed_matrix[aa_in_seed_matrix]
    nonseed_aa_matrix = nonseed_matrix[aa_in_non_seed_matrix]

    # Genome-level auxotrophy/prototrophy counts
    auxotrophy_per_genome = seed_aa_matrix.sum(axis=1)
    prototrophy_per_genome = nonseed_aa_matrix.sum(axis=1)

    summary_df = pd.DataFrame({
        'Auxotrophies': auxotrophy_per_genome,
        'Prototrophies': prototrophy_per_genome
    })

    # Counts per amino acid
    seed_counts = seed_aa_matrix.sum(axis=0)
    nonseed_counts = nonseed_aa_matrix.sum(axis=0)

    counts_df = pd.DataFrame({
        'Auxotrophy (Seed)': seed_counts,
        'Prototrophy (Non-Seed)': nonseed_counts
    })
    labels = aa_df.set_index('SEED_ID')['compound_name'].reindex(counts_df.index)
    counts_df.index = labels

    # Merge with taxonomy
    seed_aa_matrix.index = seed_aa_matrix.index.astype(str)
    metadata_df['patric_id'] = metadata_df['patric_id'].astype(str)
    merged = seed_aa_matrix.merge(
        metadata_df[['patric_id', 'gtdb_taxonomy']],
        left_index=True, right_on='patric_id', how='left'
    )
    taxonomy_df = split_and_clean_taxonomy(merged, 'gtdb_taxonomy')
    merged = pd.concat([merged, taxonomy_df], axis=1)

    return summary_df, counts_df, merged, aa_in_seed_matrix


# ------------------------
# Taxonomic aggregation
# ------------------------
def aggregate_by_taxon(merged, aa_in_seed_matrix, ranks=("phylum", "family", "genus")):
    """Compute auxotrophy counts and percentages per amino acid per taxonomic rank."""
    results = {}
    for rank in ranks:
        rank_summary = merged.groupby(rank)[aa_in_seed_matrix].sum()
        rank_counts = merged[rank].value_counts()
        rank_percentage = rank_summary.div(rank_counts, axis=0) * 100

        results[rank] = {"counts": rank_summary, "percentages": rank_percentage}
    return results


# ------------------------
# Visualization
# ------------------------
def plot_amino_acid_counts(counts_df, log_scale=False, save_path=None):
    """Barplot of auxotrophy vs prototrophy frequency per amino acid."""
    ax = counts_df.plot(kind='bar', figsize=(9, 5))
    plt.title('Auxotrophy/Prototrophy Frequency per Amino Acid')
    plt.ylabel('Number of Genomes')
    plt.xlabel('Amino Acid')

    if log_scale:
        ax.set_yscale("log")

    plt.xticks(rotation=45, ha='right')

    ax.legend(
        loc='upper center',
        bbox_to_anchor=(0.2, -0.5),   
        ncol=2,                        
        fontsize=9
    )

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Saved plot: {save_path}")
    plt.show()


def plot_summary_hist(summary_df, n_amino_acids, log_scale=False, save_path=None):
    """Histogram of auxotrophies and prototrophies per genome."""
    ax = summary_df.plot(kind='hist', bins=range(0, n_amino_acids + 2), alpha=0.7)
    plt.xlabel('Number of Amino Acids')
    plt.ylabel('Number of Genomes')

    if log_scale:
        ax.set_yscale("log")

    plt.title('Distribution of Auxotrophies and Prototrophies per Genome')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Saved plot: {save_path}")
    plt.show()


def plot_stacked(rank_summary, rank, aa_df, percentages=False, log_scale=False, top_n=20, save_path=None):
    """Stacked barplot of auxotrophy per amino acid for the top-N taxa."""
    # Compute genome counts per taxon
    genome_counts = rank_summary.sum(axis=1).sort_values(ascending=False)

    # Select top-N taxa
    top_taxa = genome_counts.head(top_n).index
    rank_summary = rank_summary.loc[top_taxa]

    # Replace SEED_IDs with compound names (and KEGG IDs)
    labels = aa_df.set_index('SEED_ID')[['compound_name']].reindex(rank_summary.columns)
    rank_summary.columns = labels.apply(lambda x: f"{x['compound_name']})", axis=1)

    # Plot
    df = rank_summary.T
    ax = df.plot(kind='bar', stacked=True, figsize=(12, 8))

    ylabel = "Percentage of Genomes" if percentages else "Number of Genomes"
    ax.set_ylabel(ylabel)
    ax.set_title(
        f"Auxotrophy (Seed) {'Percentage' if percentages else 'Counts'} "
        f"per Amino Acid per Top {top_n} {rank.capitalize()}"
    )

    if log_scale:
        ax.set_yscale("log")
    
    plt.xticks(rotation=45, ha='right')

    # Move legend below the plot
    ax.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.5),
        ncol=6,
        fontsize=8,
        title=rank.capitalize()
    )

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Saved plot: {save_path}")
    plt.show()


# ------------------------
# Main
# ------------------------
def main():
    # Load data
    seeds = load_data(SEEDS_PICKLE, filetype='pickle')
    non_seeds = load_data(NON_SEEDS_PICKLE, filetype='pickle')
    aa_df = load_data(AMINOACIDS, filetype='tsv')
    metadata_df = load_data(COMPACT_METADATA_WITH_BIOME, filetype='tsv')

    # Directory to save the plots
    plot_dir = OUTPUT_DIR / "plots" / "aminoacids"
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Analyze auxotrophy
    summary_df, counts_df, merged, aa_in_seed_matrix = analyze_auxotrophy(seeds, non_seeds, aa_df, metadata_df)

    # Taxonomic aggregation
    results = aggregate_by_taxon(merged, aa_in_seed_matrix)

    # Save outputs
    counts_df.to_csv(OUTPUT_DIR / "amino_acid_auxotrophy_prototrophy_counts.tsv", sep="\t")
    summary_df.to_csv(OUTPUT_DIR / "amino_acid_auxo_proto_per_genome.tsv", sep="\t")
    for rank, res in results.items():
        res["counts"].to_csv(OUTPUT_DIR / f"auxotrophy_counts_per_{rank}.tsv", sep="\t")
        res["percentages"].to_csv(OUTPUT_DIR / f"auxotrophy_percentages_per_{rank}.tsv", sep="\t")

    # Plots
    plot_amino_acid_counts(counts_df, log_scale=True, save_path=plot_dir/"aa_counts.png")
    plot_summary_hist(summary_df, n_amino_acids=len(aa_df), log_scale=True, save_path=plot_dir/"aa_summary_hist.png")
    for rank, res in results.items():
        plot_stacked(res["counts"], rank, aa_df, percentages=False, log_scale=True, save_path=plot_dir/f"aa_auxotrophy_counts_{rank}.png")
        plot_stacked(res["percentages"], rank, aa_df, percentages=True, log_scale=True, save_path=plot_dir/f"aa_auxotrophy_percentages_{rank}.png")


if __name__ == "__main__":
    main()
