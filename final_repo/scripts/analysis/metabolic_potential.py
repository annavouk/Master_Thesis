"""
Analyze the metabolic potential of each genome based on binary presence of seeds and non-seeds.
Calculates total seeds, total non-seeds, their ratio per genome and appends genome length from PATRIC metadata.
"""

import pandas as pd


def seeds_per_genome(df):
    """
    Count the number of seeds (value == 1) per genome.

    Args:
        df (pd.DataFrame): Binary dataframe for seeds.

    Returns:
        tuple: Updated dataframe, ID with most seeds, its count,
               ID with least seeds, its count.
    """
    df['Total_Seeds'] = (df.values == 1).sum(axis=1)

    most = df['Total_Seeds'].idxmax()
    least = df['Total_Seeds'].idxmin()

    return df, most, df.loc[most, 'Total_Seeds'], least, df.loc[least, 'Total_Seeds']


def non_seeds_per_genome(df):
    """
    Count the number of non-seeds (value == 1) per genome.

    Args:
        df (pd.DataFrame): Binary dataframe for non-seeds.

    Returns:
        tuple: Updated dataframe, ID with most non-seeds, its count,
               ID with least non-seeds, its count.
    """
    df['Total_non_Seeds'] = (df.values == 1).sum(axis=1)

    most = df['Total_non_Seeds'].idxmax()
    least = df['Total_non_Seeds'].idxmin()

    return df, most, df.loc[most, 'Total_non_Seeds'], least, df.loc[least, 'Total_non_Seeds']


if __name__ == "__main__":
    # Load binary presence/absence matrices
    seeds_df = pd.read_pickle("~/master_thesis/final_repo/input/seeds_binary_per_patric.pckl")
    non_seeds_df = pd.read_pickle("~/master_thesis/final_repo/input/non_seeds_binary_per_patric.pckl")

    # Analyze SEEDs
    seeds_df, most_seed, max_seed, least_seed, min_seed = seeds_per_genome(seeds_df)
    print(f"Most seeds: {most_seed} ({max_seed})")
    print(f"Least seeds: {least_seed} ({min_seed})")

    # Analyze non-SEEDs
    non_seeds_df, most_non, max_non, least_non, min_non = non_seeds_per_genome(non_seeds_df)
    print(f"Most non-seeds: {most_non} ({max_non})")
    print(f"Least non-seeds: {least_non} ({min_non})")

    # Keep only totals
    seeds_df = seeds_df[['Total_Seeds']].copy()
    non_seeds_df = non_seeds_df[['Total_non_Seeds']].copy()

    # Ensure PATRIC ID is index and then reset
    seeds_df.index.name = 'PATRIC'
    non_seeds_df.index.name = 'PATRIC'

    seeds_df.reset_index(inplace=True)
    non_seeds_df.reset_index(inplace=True)

    # Merge totals
    merged_df = pd.merge(seeds_df, non_seeds_df, on="PATRIC")
    merged_df['Ratio'] = merged_df['Total_Seeds'] / merged_df['Total_non_Seeds'].replace(0, pd.NA)

    # Load metadata and ensure correct type (string) for merging
    metadata_df = pd.read_csv("~/master_thesis/final_repo/metadata/compact_metadata.csv", dtype={'patric_id': str})
    metadata_df['patric_id'] = metadata_df['patric_id'].astype(str)
    merged_df['PATRIC'] = merged_df['PATRIC'].astype(str)

    metadata_df.rename(columns={'patric_id': 'PATRIC'}, inplace=True)

    # Merge genome length into final dataframe
    merged_df = pd.merge(merged_df, metadata_df[['PATRIC', 'genome_length', 'gtdb_taxonomy']], on='PATRIC', how='left')

    # Output preview
    print(merged_df.head())

    # Print stats
    print(f"Minimum ratio: {merged_df['Ratio'].min()}")
    print(f"Maximum ratio: {merged_df['Ratio'].max()}")

    # Save the result
    merged_df.to_csv("metabolic_potential_summary.csv", index=False)
