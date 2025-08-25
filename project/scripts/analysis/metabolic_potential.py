"""
Metabolic Potential Quantification (Approach 1)
Analyze the metabolic potential of each genome based on binary presence of seeds and non-seeds.
Calculate total seeds, total non-seeds, their ratio per genome and appends genome length from PATRIC metadata.
Results are saved in .csv file for visualization.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
from config import SEEDS_PICKLE, NON_SEEDS_PICKLE, COMPACT_METADATA, OUTPUT_DIR
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
    """
    Combine Total_Seeds and Total_non_Seeds counts and calculate its ratio per genome. Return final summary df.
    Output columns: patric_id, Total_Seeds, Total_non_Seeds, Ratio
    """
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

    return summary[["patric_id", "Total_Seeds", "Total_non_Seeds", "Ratio"]]


def merge_with_metadata(summary, metadata_df):
    """
    Merge summary dataframe with genome metadata for genome_length and gtdb_taxonomy.
    Output columns: patric_id, Total_Seeds, Total_non_Seeds, Ratio, genome_length, gtdb_taxonomy
    """
    # Ensure correct column name for merging
    metadata_df["patric_id"] = metadata_df["patric_id"].astype(str).str.strip()
    summary["patric_id"] = summary["patric_id"].astype(str).str.strip()

    # Merge the relevant metadata fields
    return pd.merge(
        summary,
        metadata_df[
            [
                "patric_id",
                "accession",
                "genome_length",
                "gtdb_taxonomy",
                "ncbi_taxonomy",
                "ncbi_taxid",
            ]
        ],
        on="patric_id",
        how="left",
    )


def normalize_per_mbp(df):
    """
    Add normalization per genome size (Mbp) for Seeds and non-Seeds,
    and include genome_length in Mbp for clarity.
    """
    df["genome_length"] = pd.to_numeric(df["genome_length"], errors="coerce")
    df = df[df["genome_length"] > 0].copy()
    df["genome_length_Mbp"] = df["genome_length"] / 1e6
    df["Seeds_per_Mbp"] = df["Total_Seeds"] / df["genome_length_Mbp"]
    df["Non_Seeds_per_Mbp"] = df["Total_non_Seeds"] / df["genome_length_Mbp"]
    df["Ratio_per_Mbp"] = (df["Total_Seeds"] / df["Total_non_Seeds"]) / df[
        "genome_length_Mbp"
    ]
    return df


# ------------------------
# Main
# ------------------------
def main():
    # Load
    seeds_df = load_data(SEEDS_PICKLE, filetype="pickle")
    non_seeds_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")

    metadata_df = pd.read_csv(
        COMPACT_METADATA, dtype={"patric_id": str}, low_memory=False
    )

    # Analyze SEEDs
    seeds_df, most_seed, max_seed, least_seed, min_seed = total_per_genome(seeds_df)
    print(f"Most seeds: {most_seed} ({max_seed})")
    print(f"Least seeds: {least_seed} ({min_seed})")

    # Analyze non-SEEDs
    non_seeds_df, most_non, max_non, least_non, min_non = total_per_genome(
        non_seeds_df, column="Total_non_Seeds"
    )
    print(f"Most non-seeds: {most_non} ({max_non})")
    print(f"Least non-seeds: {least_non} ({min_non})")

    # Calculate Ratio and merge with metadata
    summary = make_metabolic_potential_summary(seeds_df, non_seeds_df)
    summary_full = merge_with_metadata(summary, metadata_df)

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
    f"Ratio_per_Mbp min: {summary_norm['Ratio_per_Mbp'].min()}, max: {summary_norm['Ratio_per_Mbp'].max()}"
    )
    print(
        f"Seeds_per_Mbp min: {summary_norm['Seeds_per_Mbp'].min()}, max: {summary_norm['Seeds_per_Mbp'].max()}"
    )
    print(
        f"Non_Seeds_per_Mbp min: {summary_norm['Non_Seeds_per_Mbp'].min()}, max: {summary_norm['Non_Seeds_per_Mbp'].max()}"
    )

    # Save output (normalized with per Mbp columns)
    out_path = OUTPUT_DIR / "metabolic_potential_summary_per_mbp.csv"
    summary_norm.to_csv(out_path, index=False)
    print(f"Saved metabolic potential summary (with per Mbp columns) to: {out_path}")


if __name__ == "__main__":
    main()
