"""
Metabolic Potential Quantification (Approach 1)
Analyze the metabolic potential of each genome based on binary presence of seeds and non-seeds.
Calculates total seeds, total non-seeds, their ratio per genome and appends genome length from PATRIC metadata.
Results are saved in .csv file for visualization.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
from config import SEEDS_PICKLE, NON_SEEDS_PICKLE, COMPACT_METADATA, OUTPUT_DIR
from utils import load_data


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


def make_metabolic_potential_summary(seeds_df, non_seeds_df):
    """
    Combine and calculate ratio. Return final summary df.
    Output columns: patric_id, Total_Seeds, Total_non_Seeds, Ratio
    """
    # Keep Total number of Seeds and non-Seeds
    seeds = seeds_df[['Total_Seeds']]
    nonseeds = non_seeds_df[['Total_non_Seeds']]

    # Ensure PATRIC_IDs is index and then reset as column
    seeds.index.name = 'PATRIC'
    nonseeds.index.name = 'PATRIC'
    seeds = seeds.reset_index()
    nonseeds = nonseeds.reset_index()

    # Merge
    summary = pd.merge(seeds, nonseeds, on="PATRIC")
    summary['Ratio'] = summary['Total_Seeds'] / summary['Total_non_Seeds'].replace(0, pd.NA)

    summary = summary.rename(columns={'PATRIC': 'patric_id'})

    return summary[['patric_id', 'Total_Seeds', 'Total_non_Seeds', 'Ratio']]


def merge_with_metadata(summary, metadata_df):
    """
    Merge summary dataframe with genome metadata for genome_length and gtdb_taxonomy.

    Args:
        summary (pd.DataFrame): DataFrame with at least 'patric_id' column.
        metadata_df (pd.DataFrame): Genome metadata (must have 'patric_id' or 'PATRIC', 'genome_length', 'gtdb_taxonomy').

    Returns:
        pd.DataFrame: merged DataFrame with extra columns for genome length and taxonomy.
    """
    # Ensure correct column name for merging
    metadata_df['patric_id'] = metadata_df['patric_id'].astype(str).str.strip()
    summary['patric_id'] = summary['patric_id'].astype(str).str.strip()
    print(repr(metadata_df[metadata_df['patric_id'].str.startswith('100233')]['patric_id'].tolist()))

    # Merge the relevant metadata fields
    return pd.merge(
        summary,
        metadata_df[["patric_id", "genome_length", "gtdb_taxonomy"]],
        on="patric_id",
        how="left"
    )

def main():
    # Load
    seeds_df = load_data(SEEDS_PICKLE, filetype="pickle")
    non_seeds_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")
    #metadata_df = load_data(COMPACT_METADATA, filetype="csv")
    metadata_df = pd.read_csv(COMPACT_METADATA, dtype={'patric_id': str})

    # Analyze SEEDs
    seeds_df, most_seed, max_seed, least_seed, min_seed = seeds_per_genome(seeds_df)
    print(f"Most seeds: {most_seed} ({max_seed})")
    print(f"Least seeds: {least_seed} ({min_seed})")

    # Analyze non-SEEDs
    non_seeds_df, most_non, max_non, least_non, min_non = non_seeds_per_genome(non_seeds_df)
    print(f"Most non-seeds: {most_non} ({max_non})")
    print(f"Least non-seeds: {least_non} ({min_non})")

    # Calculate Ratio and merge with metadata
    summary =  make_metabolic_potential_summary(seeds_df, non_seeds_df)
    summary.to_csv(OUTPUT_DIR / "metabolic_potential_summary_no_metadata.csv", index=False)
    summary_full = merge_with_metadata(summary, metadata_df)

    # Print summary stats
    print(summary_full.head())
    print(f"Ratio min: {summary_full['Ratio'].min()}, max: {summary_full['Ratio'].max()}")

    # Save output
    out_path = OUTPUT_DIR / "metabolic_potential_summary.csv"
    summary_full.to_csv(out_path, index=False)
    print(f"Saved metabolic potential summary to: {out_path}")


if __name__ == "__main__":
    main()
