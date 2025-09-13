"""
Merge genome metadata from PATRIC and GTDB.

- Cleans and harmonizes accession IDs
- Performs dual merge (NCBI vs GTDB accession)
- Adds environmental annotations from PREGO (lit_envs, meta_envs, all_envs)
- Outputs merged table with taxonomic and genome information
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd

from config import (
    PATRIC_METADATA_JSON,  # input (previously fetched)
    PREGO_ENVIRONMENTS_JSON,  # input (previously fetched)
    GTDB_METADATA_TSV,  # input (downloaded and merged)
    PREGO_INPUT_CSV,  # output (after dual merge, to be used for PREGO)
    COMPACT_METADATA_TSV,  # output
)
from utils import load_data


# ------------------------
# Cleaning helpers
# ------------------------
def clean_patric_accessions(patric_df):
    """Clean and add root accession (no version) for merge."""
    patric_df_filtered = patric_df[["patric_id", "assembly_accession"]].copy()
    patric_df_filtered["root_accession"] = (
        patric_df_filtered["assembly_accession"]
        .fillna("")
        .str.replace(r"^(GCA_|GCF_)", "", regex=True)
        .str.split(".")
        .str[0]
        .astype(str)
    )
    return patric_df_filtered


def clean_gtdb_metadata(gtdb_df):
    """Clean GTDB: create root accessions from both NCBI accession and GTDB accession."""
    columns_to_keep = [
        "accession",
        "genome_size",
        "gtdb_taxonomy",
        "ncbi_genbank_assembly_accession",
        "ncbi_taxid",
        "checkm_completeness",
        "checkm_contamination",
        "gtdb_representative",
    ]
    gtdb_df_filtered = gtdb_df[columns_to_keep].copy()

    # From NCBI
    gtdb_df_filtered["root_ncbi_accession"] = (
        gtdb_df_filtered["ncbi_genbank_assembly_accession"]
        .fillna("")
        .str.replace(r"^(GCA_|GCF_)", "", regex=True)
        .str.split(".")
        .str[0]
        .astype(str)
    )

    # From GTDB 'accession'
    gtdb_df_filtered["root_gtdb_accession"] = (
        gtdb_df_filtered["accession"]
        .fillna("")
        .str.replace(r"^(RS_GCF_|GB_GCA_|GB_GCF_|GCA_|GCF_)", "", regex=True)
        .str.split(".")
        .str[0]
        .astype(str)
    )

    return gtdb_df_filtered


# ------------------------
# PREGO enrichment
# ------------------------
def add_prego_envs(df, prego):
    """Add lit_envs, meta_envs, all_envs columns from PREGO to metadata table."""
    taxid_str = df["ncbi_taxid"].astype(str)

    # Fetch lit and meta environments as lists
    lit_envs = taxid_str.map(
        lambda taxid: prego.get(taxid, {}).get("lit_envs", []) or []
    )
    meta_envs = taxid_str.map(
        lambda taxid: prego.get(taxid, {}).get("meta_envs", []) or []
    )

    # Convert lists to string columns
    df["lit_envs"] = lit_envs.map(lambda x: "; ".join(x))
    df["meta_envs"] = meta_envs.map(lambda x: "; ".join(x))

    # Union of lit and meta envs
    df["all_envs"] = [
        "; ".join(sorted(set(lit) | set(meta)))
        for lit, meta in zip(lit_envs, meta_envs)
    ]

    return df


# ------------------------
# Merging logic
# ------------------------
def dual_merge(patric_df, gtdb_df, output_path=None):
    """Merge first on root_accession (PATRIC) vs root_ncbi_accession (GTDB), then unmatched on root_accession vs root_gtdb_accession (GTDB).
    Returns single, deduplicated merged DataFrame."""
    # First merge: NCBI accession
    merged_1 = pd.merge(
        patric_df,
        gtdb_df,
        left_on="root_accession",
        right_on="root_ncbi_accession",
        how="left",
        suffixes=("_patric", "_gtdb"),
    )
    unmatched = merged_1[merged_1["gtdb_taxonomy"].isna()]
    print(
        f"After 1st merge (NCBI): unmatched genomes: {unmatched.shape[0]} / {patric_df.shape[0]}"
    )

    # Only unmatched go for 2nd merge (GTDB accession)
    unmatched_patric = unmatched[["patric_id", "assembly_accession", "root_accession"]]
    merged_2 = pd.merge(
        unmatched_patric,
        gtdb_df,
        left_on="root_accession",
        right_on="root_gtdb_accession",
        how="left",
    )
    print(
        f"After 2nd merge (GTDB accession): rescued: {merged_2[~merged_2['gtdb_taxonomy'].isna()].shape[0]}"
    )

    # Combine results: matched from 1st, rescued from 2nd, remove duplicates
    matched_1 = merged_1[~merged_1["gtdb_taxonomy"].isna()]
    matched_2 = merged_2[~merged_2["gtdb_taxonomy"].isna()]
    merged_final = pd.concat([matched_1, matched_2], ignore_index=True)
    merged_final = merged_final.drop_duplicates(subset=["patric_id"])
    print(
        f"Total matched after dual merge: {merged_final.shape[0]} / {patric_df.shape[0]}"
    )

    # Save
    if output_path:
        merged_final.to_csv(output_path, index=False)
        print(f"Merged metadata saved to: {output_path}")

    # Unmatched report
    unmatched_final = set(patric_df["patric_id"]) - set(merged_final["patric_id"])
    if unmatched_final:
        print("Remaining unmatched PATRIC assemblies:")
        print(list(unmatched_final)[:10])
    else:
        print("All PATRIC assemblies matched after dual merge!")

    return merged_final


# ------------------------
# Helper
# ------------------------
def genome_size_per_mbp(df):
    """Convert genome size (bp) -> (Mbp)."""
    df = df.rename(columns={"genome_size": "genome_size_bp"})
    df["genome_size_bp"] = pd.to_numeric(df["genome_size_bp"], errors="coerce")
    cdf = df[df["genome_size_bp"] > 0].copy()
    cdf["genome_size_Mbp"] = cdf["genome_size_bp"] / 1e6

    return cdf


# ------------------------
# Main
# ------------------------
def main():
    # Load data
    patric_data = load_data(PATRIC_METADATA_JSON, filetype="json")
    patric_df = pd.DataFrame.from_dict(patric_data, orient="index").reset_index()
    patric_df.rename(columns={"index": "patric_id"}, inplace=True)
    gtdb_df = load_data(GTDB_METADATA_TSV, filetype="tsv")
    prego = load_data(PREGO_ENVIRONMENTS_JSON, filetype="json")

    # Create root accession for merge and keep only columns of interest
    patric_clean = clean_patric_accessions(patric_df)
    gtdb_clean = clean_gtdb_metadata(gtdb_df)

    # Merge
    merged = dual_merge(patric_clean, gtdb_clean)
    cols_to_keep = [
        "patric_id",
        "genome_size",
        "gtdb_taxonomy",
        "ncbi_taxid",
        "accession",
        "assembly_accession",
        "ncbi_genbank_assembly_accession",
        "checkm_completeness",
        "checkm_contamination",
        "gtdb_representative",
    ]
    merged = merged[cols_to_keep]
    merged.to_csv(PREGO_INPUT_CSV, index=False)  # To fetch PREGO envs

    print(f"Merged rows:  {len(merged):,}")
    missing_taxid = merged["ncbi_taxid"].isna().sum()
    print(f"Rows without ncbi_taxid (so no PREGO): {missing_taxid:,}")

    # Add PREGO environments
    merged = add_prego_envs(merged, prego)

    nonempty_all = (merged["all_envs"].str.strip() != "").sum()
    print(
        f"all_envs non-empty: {nonempty_all}/{len(merged)} ({nonempty_all/len(merged):.1%})"
    )

    # bp -> Mbp
    merged = genome_size_per_mbp(merged)

    # Save final
    merged.to_csv(COMPACT_METADATA_TSV, sep="\t", index=False)
    print(f"Final merged metadata + environments saved to: {COMPACT_METADATA_TSV}")


if __name__ == "__main__":
    main()
