import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
from config import PATRIC_METADATA_JSON, GTDB_METADATA, COMPACT_METADATA
from utils import load_data

def clean_patric_accessions(df):
    """
    Clean and add root accession (no version) for merge.
    """
    df_filtered = df[['patric_id', 'assembly_accession', 'genome_length']].copy()
    df_filtered['assembly_accession'] = df_filtered['assembly_accession'].str.replace(
        r'^(GCA_|GCF_)', '', regex=True
    )
    df_filtered['root_assembly'] = df_filtered['assembly_accession'].str.split('.').str[0]
    return df_filtered

def clean_gtdb_metadata(df):
    """
    Clean GTDB: create root accessions from both NCBI accession and GTDB accession.
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
    # From NCBI
    df_filtered['assembly_accession_ncbi'] = df_filtered['ncbi_genbank_assembly_accession'].astype(str).str.replace(
        r'^(GCA_|GCF_)', '', regex=True
    )
    df_filtered['root_assembly_ncbi'] = df_filtered['assembly_accession_ncbi'].str.split('.').str[0]
    # From GTDB 'accession'
    df_filtered['assembly_accession_acc'] = df_filtered['accession'].astype(str).str.replace(
        r'^(RS_GCF_|GB_GCA_|GB_GCF_|GCA_|GCF_)', '', regex=True
    )
    df_filtered['root_assembly_acc'] = df_filtered['assembly_accession_acc'].str.split('.').str[0]
    return df_filtered

def dual_merge(patric_df, gtdb_df, output_path=None):
    """
    Merge first on root_assembly vs root_assembly_ncbi,
    then unmatched on root_assembly vs root_assembly_acc.
    Returns single, deduplicated merged DataFrame.
    """
    # First merge: PATRIC vs GTDB by NCBI accession
    merged_1 = pd.merge(
        patric_df, gtdb_df,
        left_on='root_assembly', right_on='root_assembly_ncbi',
        how='left', suffixes=('_patric', '_gtdb')
    )
    unmatched = merged_1[merged_1['gtdb_taxonomy'].isna()]
    print(f"After 1st merge (NCBI): unmatched genomes: {unmatched.shape[0]} / {patric_df.shape[0]}")

    # Only unmatched go for 2nd merge (using GTDB accession)
    unmatched_patric = unmatched[['patric_id', 'assembly_accession', 'genome_length', 'root_assembly']]
    merged_2 = pd.merge(
        unmatched_patric, gtdb_df,
        left_on='root_assembly', right_on='root_assembly_acc',
        how='left'
    )
    print(f"After 2nd merge (GTDB accession): rescued: {merged_2[~merged_2['gtdb_taxonomy'].isna()].shape[0]}")

    # Combine results: matched from 1st, rescued from 2nd, remove duplicates
    matched_1 = merged_1[~merged_1['gtdb_taxonomy'].isna()]
    matched_2 = merged_2[~merged_2['gtdb_taxonomy'].isna()]
    merged_final = pd.concat([matched_1, matched_2], ignore_index=True)
    merged_final = merged_final.drop_duplicates(subset=['patric_id'])
    print(f"Total matched after dual merge: {merged_final.shape[0]} / {patric_df.shape[0]}")

    # Save
    if output_path:
        merged_final.to_csv(output_path, index=False)
        print(f"Merged metadata saved to: {output_path}")

    # Unmatched report
    unmatched_final = set(patric_df['patric_id']) - set(merged_final['patric_id'])
    if unmatched_final:
        print("Remaining unmatched PATRIC assemblies:")
        print(list(unmatched_final)[:10])
    else:
        print("All PATRIC assemblies matched after dual merge!")

    return merged_final

def main():
    patric_data = load_data(PATRIC_METADATA_JSON, filetype='json')
    patric_df = pd.DataFrame.from_dict(patric_data, orient='index').reset_index()
    patric_df.rename(columns={'index': 'patric_id'}, inplace=True)
    gtdb_df = load_data(GTDB_METADATA, filetype='tsv')

    patric_clean = clean_patric_accessions(patric_df)
    gtdb_clean = clean_gtdb_metadata(gtdb_df)

    dual_merge(patric_clean, gtdb_clean, output_path=COMPACT_METADATA)

if __name__ == "__main__":
    main()
