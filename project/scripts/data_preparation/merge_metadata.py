"""
Merge PATRIC metadata with GTDB taxonomy information based on assembly accessions.
Creates a compact metadata table for downstream analysis.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import json
from config import PATRIC_METADATA_JSON, GTDB_METADATA, COMPACT_METADATA
from utils import load_data


def clean_patric_accessions(df):
    """
    Filters relevant columns and cleans assembly accession format.

    Parameters:
    - df: pandas.DataFrame with at least ['patric_id', 'assembly_accession', 'genome_length']

    Returns:
    - Cleaned DataFrame with filtered and formatted assembly accessions
    """
    df_filtered = df[['patric_id', 'assembly_accession', 'genome_length']].copy()
    df_filtered['assembly_accession'] = df_filtered['assembly_accession'].str.replace(
        r'^(GCA_|GCF_)', '', regex=True
    )
    return df_filtered

def clean_gtdb_metadata(df):
    """
    Filters and cleans GTDB metadata by normalizing assembly accessions.

    Parameters:
    - df: pandas.DataFrame with GTDB metadata

    Returns:
    - Cleaned DataFrame with standardized 'assembly_accession' column
    """
    columns_to_keep = [
        'accession',
        'genome_size',
        'gtdb_taxonomy',
        'ncbi_genbank_assembly_accession',
        'ncbi_taxonomy',
        'checkm_completeness',
        'checkm_contamination',
        'gtdb_genome_representative'
    ]
    df_filtered = df[columns_to_keep].copy()

    # Remove NCBI-style prefixes if present
    df_filtered['ncbi_genbank_assembly_accession'] = df_filtered['ncbi_genbank_assembly_accession'].str.replace(
        r'^(GCA_|GCF_)', '', regex=True)

    # Rename for merging
    df_filtered.rename(columns={'ncbi_genbank_assembly_accession': 'assembly_accession'}, inplace=True)

    return df_filtered

def merge_metadata(patric_df, gtdb_df, output_path=None):
    """
    Merge cleaned PATRIC and GTDB metadata on 'assembly_accession'.

    Parameters:
    - patric_df: DataFrame with cleaned PATRIC metadata
    - gtdb_df: DataFrame with cleaned GTDB metadata
    - output_path: Optional path to save the merged DataFrame as CSV

    Returns:
    - merged_df: The merged pandas DataFrame
    """
    merged_df = pd.merge(patric_df, gtdb_df, on='assembly_accession', how='left')

    # Report unmatched genomes
    unmatched = merged_df['gtdb_taxonomy'].isna().sum()
    print(f"Unmatched genomes: {unmatched} out of {len(merged_df)}")

    if output_path:
        merged_df.to_csv(output_path, index=False)
        print(f"Merged metadata saved to: {output_path}")

    return merged_df

def main():
    # Load raw metadata
    patric_data = load_data(PATRIC_METADATA_JSON, filetype='json')
    patric_df = pd.DataFrame.from_dict(patric_data, orient='index').reset_index()
    patric_df.rename(columns={'index': 'patric_id'}, inplace=True)

    gtdb_df = load_data(GTDB_METADATA, filetype='tsv')


    # Clean & merge
    patric_clean = clean_patric_accessions(patric_df)
    gtdb_clean = clean_gtdb_metadata(gtdb_df)

    print("PATRIC unique assembly_accession:", patric_clean['assembly_accession'].nunique())
    print("GTDB unique assembly_accession:", gtdb_clean['assembly_accession'].nunique())
    print("PATRIC duplicates:", patric_clean['assembly_accession'].duplicated().sum())
    print("GTDB duplicates:", gtdb_clean['assembly_accession'].duplicated().sum())

    merged_df = merge_metadata(patric_clean, gtdb_clean, output_path=COMPACT_METADATA)

    print("Rows in merged:", merged_df.shape[0])
    print("Unique patric_id in merged:", merged_df['patric_id'].nunique())
    print("patric_id duplicates:", merged_df['patric_id'].duplicated().sum())

    unmatched = merged_df[merged_df['gtdb_taxonomy'].isna()]
    print(unmatched[['patric_id', 'assembly_accession']].head(20))


if __name__ == "__main__":
    main()
