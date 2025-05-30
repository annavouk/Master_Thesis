"""
Analyze which compounds, come with higher seeds/non-seeds frequency, based on binary presence.
Calculates total genomes count per seed/non-seed, combines them and appends the KEGG modules the compound is present in.
"""

import pandas as pd
from collections import defaultdict


def genomes_per_seed(df):
    """
    Count the number of genomes (value == 1) per seed compound.

    Args:
        df (pd.DataFrame): Binary dataframe for seeds.

    Returns:
        tuple: Updated dataframe, seed compound with most genomes, its count,
               seed compound with least genomes, its count.
    """
    genomes = df.sum()
    most_genomes = genomes.idxmax()
    least_genomes = genomes.idxmin()
    max_count = genomes.max()
    min_count = genomes.min()

    sum_row = pd.DataFrame(genomes).T
    sum_row.index = ['Total_Genomes_per_Seed']
    updated_df = pd.concat([df, sum_row])

    return (genomes, most_genomes, max_count, least_genomes, min_count, updated_df)


def genomes_per_non_seed(df):
    """
    Count the number of genomes (value == 1) per non-seed compound.

    Args:
        df (pd.DataFrame): Binary dataframe for non-seeds.

    Returns:
        tuple: Updated dataframe, non-seed compound with most genomes, its count,
               non-seed compound with least genomes, its count.
    """
    # Count 1s per column
    genomes = df.sum()
    most_genomes = genomes.idxmax()
    least_genomes = genomes.idxmin()
    max_count = genomes.max()
    min_count = genomes.min()

    sum_row = pd.DataFrame(genomes).T
    sum_row.index = ['Total_Genomes_per_non_Seed']
    updated_df = pd.concat([df, sum_row])

    return (genomes, most_genomes, max_count, least_genomes, min_count, updated_df)


def combine_seed_and_nonseed_counts(seed_counts, non_seed_counts):
    """
    Combine seed and non-seed genome counts into a single DataFrame.

    Args:
        seed_counts (pd.Series): Genome counts per seed compound.
        non_seed_counts (pd.Series): Genome counts per non-seed compound.

    Returns:
        pd.DataFrame: Combined DataFrame with both counts, filling missing values with 0.
    """
    seed_counts.name = "seed_genome_count"
    non_seed_counts.name = "non_seed_genome_count"

    combined_df = pd.concat([seed_counts, non_seed_counts], axis=1).fillna(0).astype(int)
    return combined_df


if __name__ == "__main__":
    # Load binary presence/absence matrices
    seeds_df = pd.read_pickle("/home/annavouk/master_thesis/final_repo/input/seeds_binary_per_patric.pckl")
    non_seeds_df = pd.read_pickle("/home/annavouk/master_thesis/final_repo/input/non_seeds_binary_per_patric.pckl")

    # Analyze seeds
    (
        seed_counts,
        most_seed,
        max_seed,
        least_seed,
        min_seed,
        updated_seeds_df
    ) = genomes_per_seed(seeds_df)

    # Analyze non-seeds
    (
        non_seed_counts,
        most_non_seed,
        max_non_seed,
        least_non_seed,
        min_non_seed,
        updated_non_seeds_df
    ) = genomes_per_non_seed(non_seeds_df)

    # Combine results
    combined_counts_df = combine_seed_and_nonseed_counts(seed_counts, non_seed_counts)
    combined_counts_df.reset_index(inplace=True)
    combined_counts_df.rename(columns={"index": "ModelSEED ID"}, inplace=True)

    # Load KEGG mapping
    column_names = ["ModelSEED ID", "KEGG ID", "KEGG Module"]
    kegg_map = pd.read_csv("~/master_thesis/final_repo/seedId_keggId_module.tsv", sep='\t', names=column_names, low_memory=False)

    # Group KEGG modules per ModelSEED ID
    module_info = kegg_map.groupby('ModelSEED ID').agg({
    'KEGG ID': 'first',
    'KEGG Module': lambda x: sorted(set(x))
    }).reset_index()

    # Merge back into main table
    df_merged = pd.merge(combined_counts_df.drop(columns=['KEGG ID', 'KEGG Module'], errors='ignore'), module_info, on='ModelSEED ID', how='left')

    # Save
    df_merged.to_csv('compounds_seed_nonseed_with_kegg.csv', index=False)

    """
    This script calculates metabolic complementarity between genomes based on their seed (required) and non-seed (produced) metabolites.
    Given two binary DataFrames, it outputs a matrix where each cell [i, j] shows how many metabolites genome i can provide to genome j,
    helping identify potential metabolic givers and receivers.
    """
    # Create dict to store output
#    nonseed_to_seed_stats = {}

    # Loop through each genome
#    for genome in non_seeds_df.index:
        # Get this genome's non-seed compounds
#        genome_nonseeds = non_seeds_df.loc[genome]
#        genome_nonseeds = genome_nonseeds[genome_nonseeds == 1].index.tolist()

        # Filter compounds to only those present in seeds_df
#        genome_nonseeds = [cpd for cpd in genome_nonseeds if cpd in seeds_df.columns]

#        if not genome_nonseeds:
#            continue  # skip genomes with no valid compounds

        # Get the subset of the seed matrix for those compounds
#        seeds_subset = seeds_df[genome_nonseeds]

        # Remove the current genome to avoid self-comparison
#        seeds_subset = seeds_subset.drop(index=genome, errors='ignore')

        # Find overlaps
#        compound_to_genomes = {}
#        for compound in seeds_subset.columns:
#            overlapping_genomes = seeds_subset.index[seeds_subset[compound] == 1].tolist()
#            if overlapping_genomes:
#                compound_to_genomes[compound] = {
#                    "count": len(overlapping_genomes),
#                    "genomes": overlapping_genomes
#                }

#        nonseed_to_seed_stats[genome] = compound_to_genomes

#    records = []

#    for genome, compounds in nonseed_to_seed_stats.items():
#        for compound, info in compounds.items():
#            records.append({
#                "genome": genome,
#                "compound": compound,
#                "num_other_genomes_seeding_it": info["count"],
#                "seeder_genomes": ";".join(info["genomes"])
#            })

#    nonseed_seed_df = pd.DataFrame(records)
#    print(non_seed_df)
