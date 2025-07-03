"""
Merge PATRIC metadata with GTDB taxonomy information based on assembly accessions.
Creates a compact metadata table for downstream analysis.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
from config import PATRIC_METADATA_JSON, GTDB_METADATA, COMPACT_METADATA
from utils import load_data

def clean_patric_accessions(df):
    """
    Filters relevant columns, cleans assembly accession format,
    and adds root accession (without version) for merge.
    """
    df_filtered = df[['patric_id', 'assembly_accession', 'genome_length']].copy()
    # Remove prefixes GCA_/GCF_
    df_filtered['assembly_accession'] = df_filtered['assembly_accession'].str.replace(
        r'^(GCA_|GCF_)', '', regex=True
    )
    # Get root accession (remove version)
    df_filtered['root_assembly'] = df_filtered['assembly_accession'].str.split('.').str[0]
    return df_filtered

def clean_gtdb_metadata(df):
    """
    Filters and cleans GTDB metadata, normalizes assembly accessions,
    and adds root accession (without version) for merge.
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
    # Remove prefix from NCBI accession
    df_filtered['assembly_accession'] = df_filtered['ncbi_genbank_assembly_accession'].astype(str).str.replace(
        r'^(GCA_|GCF_)', '', regex=True
    )
    # Get root accession (remove version)
    df_filtered['root_assembly'] = df_filtered['assembly_accession'].str.split('.').str[0]
    return df_filtered

def merge_metadata(patric_df, gtdb_df, output_path=None):
    """
    Merge on 'root_assembly', keeping all PATRIC rows.
    """
    merged_df = pd.merge(
        patric_df, gtdb_df,
        on='root_assembly',
        how='left',
        suffixes=('_patric', '_gtdb')
    )
    unmatched = merged_df['gtdb_taxonomy'].isna().sum()
    print(f"Unmatched genomes: {unmatched} out of {len(merged_df)}")

    if output_path:
        merged_df.to_csv(output_path, index=False)
        print(f"Merged metadata saved to: {output_path}")

    return merged_df

def main():
    # Load and clean metadata
    patric_data = load_data(PATRIC_METADATA_JSON, filetype='json')
    patric_df = pd.DataFrame.from_dict(patric_data, orient='index').reset_index()
    patric_df.rename(columns={'index': 'patric_id'}, inplace=True)
    gtdb_df = load_data(GTDB_METADATA, filetype='tsv')

    patric_clean = clean_patric_accessions(patric_df)
    gtdb_clean = clean_gtdb_metadata(gtdb_df)

    print("PATRIC unique root assemblies:", patric_clean['root_assembly'].nunique())
    print("GTDB unique root assemblies:", gtdb_clean['root_assembly'].nunique())

    merged_df = merge_metadata(patric_clean, gtdb_clean, output_path=COMPACT_METADATA)

    print("Rows in merged:", merged_df.shape[0])
    print("Unique patric_id in merged:", merged_df['patric_id'].nunique())
    print("patric_id duplicates:", merged_df['patric_id'].duplicated().sum())

    unmatched = merged_df[merged_df['gtdb_taxonomy'].isna()]
    if not unmatched.empty:
        print("First 20 unmatched PATRIC assemblies:")
        print(unmatched[['patric_id', 'assembly_accession_patric']].head(20))
    else:
        print("All PATRIC assemblies matched to GTDB.")


if __name__ == "__main__":
    main()
