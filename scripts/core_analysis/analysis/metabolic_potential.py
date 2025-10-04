"""
Metabolic Potential Quantification

Analyze the metabolic potential of each genome based on binary presence of seeds and non-seeds.
Calculate total Seeds, total non-Seeds, their ratio per genome and appends genome size (Mbp) from PATRIC metadata.
Results are saved in .tsv file for visualization.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

import pandas as pd

from config import( 
    SEEDS_PICKLE,    # input
    NON_SEEDS_PICKLE,    # input
    COMPACT_METADATA_TSV,    # input
    METABOLIC_POTENTIAL_TSV,    # output
    )
from utils import load_data


# ------------------------
# Helpers
# ------------------------
def total_per_genome(df, column="Total_Seeds"):
    """Count the number of seeds or non-seeds (value == 1) per genome and add it as a column and return min/max stats."""
    df[column] = (df.values == 1).sum(axis=1)

    most = df[column].idxmax()
    least = df[column].idxmin()

    return df, most, df.loc[most, column], least, df.loc[least, column]


def make_metabolic_potential_summary(seeds_df, non_seeds_df):
    """Combine Total_Seeds and Total_non_Seeds counts and calculate its ratio per genome."""
    # Keep Total number of Seeds and non-Seeds
    seeds = seeds_df[["Total_Seeds"]]
    nonseeds = non_seeds_df[["Total_non_Seeds"]]

    # Ensure PATRIC_IDs is index and then reset as column
    seeds.index.name = "PATRIC"
    nonseeds.index.name = "PATRIC"
    seeds = seeds.reset_index()
    nonseeds = nonseeds.reset_index()

    # Merge
    summary = pd.merge(seeds, nonseeds, on="PATRIC")
    summary["Ratio"] = summary["Total_Seeds"] / summary["Total_non_Seeds"].replace(
        0, pd.NA
    )

    summary = summary.rename(columns={"PATRIC": "patric_id"})
    summary = summary[["patric_id", "Total_Seeds", "Total_non_Seeds", "Ratio"]]

    return summary


def merge_with_metadata(summary, metadata_df):
    """Merge summary with genome_size_Mbp and gtdb_taxonomy."""
    # Ensure correct column name for merging
    metadata_df["patric_id"] = metadata_df["patric_id"].astype(str).str.strip()
    summary["patric_id"] = summary["patric_id"].astype(str).str.strip()

    # Merge the relevant metadata fields
    merged = pd.merge(
        summary,
        metadata_df[["patric_id", "genome_size_Mbp", "gtdb_taxonomy"]],
        on="patric_id",
        how="left",
    )
    merged = merged[
        ["patric_id", "Total_Seeds", "Total_non_Seeds", "Ratio",
         "genome_size_Mbp", "gtdb_taxonomy"]
    ]

    return merged


def normalize_per_mbp(df):
    """Add normalization per genome size (Mbp) for Seeds and non-Seeds."""
    df["Seeds_per_Mbp"] = df["Total_Seeds"] / df["genome_size_Mbp"].replace(
        0, pd.NA
    )
    df["Non_Seeds_per_Mbp"] = df["Total_non_Seeds"] / df["genome_size_Mbp"].replace(
        0, pd.NA
    )

    return df


# ------------------------
# Main
# ------------------------
def main():
    # Load
    seeds_df = load_data(SEEDS_PICKLE, filetype="pickle")
    non_seeds_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")

    metadata_df = load_data(COMPACT_METADATA_TSV, filetype="tsv", dtype={"patric_id": str})

    # Analyze Seeds
    seeds_df, most_seed, max_seed, least_seed, min_seed = total_per_genome(seeds_df)
    print(f"Most Seeds: {most_seed} ({max_seed})")
    print(f"Least Seeds: {least_seed} ({min_seed})")

    # Analyze non-Seeds
    non_seeds_df, most_non, max_non, least_non, min_non = total_per_genome(
        non_seeds_df, column="Total_non_Seeds"
    )
    print(f"Most non-Seeds: {most_non} ({max_non})")
    print(f"Least non-Seeds: {least_non} ({min_non})")

    # Calculate Ratio and merge with metadata
    summary = make_metabolic_potential_summary(seeds_df, non_seeds_df)
    summary_full = merge_with_metadata(summary, metadata_df)
    summary_full["genome_size_Mbp"] = pd.to_numeric(
    summary_full["genome_size_Mbp"], errors="coerce"
    )


    # Print summary stats
    print(summary_full.head())
    print(
        f"Ratio min: {summary_full['Ratio'].min()}, max: {summary_full['Ratio'].max()}"
    )

    # Normalization per Mbp
    summary_norm = normalize_per_mbp(summary_full)

    # Print preview
    print(summary_norm.head())
    print(
        f"Seeds_per_Mbp min: {summary_norm['Seeds_per_Mbp'].min()}, max: {summary_norm['Seeds_per_Mbp'].max()}"
    )
    print(
        f"Non_Seeds_per_Mbp min: {summary_norm['Non_Seeds_per_Mbp'].min()}, max: {summary_norm['Non_Seeds_per_Mbp'].max()}"
    )

    # Save output (normalized with per Mbp columns)
    out_file = Path(METABOLIC_POTENTIAL_TSV)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    summary_norm.to_csv(METABOLIC_POTENTIAL_TSV, sep='\t', index=False)
    print(f"Saved metabolic potential summary to: {out_file}")


if __name__ == "__main__":
    main()
