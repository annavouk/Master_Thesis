"""
Metabolic Potential - Genome-by-Pathway Functional Coverage Analysis (Approach 2)

For each genome, calculates the coverage (% of compounds present as non-seeds) for each KEGG pathway.
"""
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import numpy as np
from config import SEEDS_PICKLE, NON_SEEDS_PICKLE, COMPACT_METADATA, COMPOUND_SUMMARY_TSV, OUTPUT_DIR
from utils import load_data
from collections import defaultdict


def get_pathway_to_compounds_map(cpd_df):
    """
    Returns a mapping: KEGG pathway ID -> set of SEED_IDs (compounds)
    participating in each pathway.

    Args:
        cpd_df (pd.DataFrame): Compound summary dataframe with columns 'SEED_ID' and 'KEGG_pathways'.

    Returns:
        dict: {KEGG_pathway: set(SEED_IDs)}
    """
    pathway_to_cpds = dict()
    for _, row in cpd_df.iterrows():
        if pd.isna(row.get('KEGG_pathways', None)):
            continue
        for pw in str(row['KEGG_pathways']).split(','):
            pw = pw.strip()
            if pw:
                pathway_to_cpds.setdefault(pw, set()).add(row['SEED_ID'])
    return pathway_to_cpds

def compute_genome_pathway_coverage(nonseed_df, pathway_to_cpds):
    """
    For each genome and each KEGG pathway, computes:
      - n_covered: Number of pathway compounds present as non-seeds in the genome
      - n_required: Total number of compounds in the pathway
      - coverage: Fraction (n_covered / n_required)

    Args:
        nonseed_df (pd.DataFrame): Binary df (genome x SEED_IDs), 1 if genome has compound as non-seed.
        pathway_to_cpds (dict): KEGG pathway ID -> set of SEED_IDs

    Returns:
        pd.DataFrame: Table with [patric_id, KEGG_pathway, n_required, n_covered, coverage]
    """
    results = []
    for genome in nonseed_df.index:
        # Non-seed compounds for this genome
        genome_cpds = set(nonseed_df.columns[nonseed_df.loc[genome] == 1])
        for pw, pw_cpds in pathway_to_cpds.items():
            n_required = len(pw_cpds)
            if n_required == 0: continue
            n_covered = len(genome_cpds & pw_cpds)
            frac = n_covered / n_required
            results.append({
                "patric_id": genome,
                "KEGG_pathway": pw,
                "n_required": n_required,
                "n_covered": n_covered,
                "coverage": frac
            })
    return pd.DataFrame(results)


def main():
    # Load data
    #seed_df = load_data(SEEDS_PICKLE, filetype="pickle")
    nonseed_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")
    summary_df = load_data(COMPOUND_SUMMARY_TSV, filetype="tsv")
    
    # Output directory
    OUTPUT = Path(OUTPUT_DIR) / "pathway_coverage_results.tsv"

    # Get pathway to compounds mapping
    pathway_to_cpds = get_pathway_to_compounds_map(summary_df)
    print(pathway_to_cpds)

    # Coverage table
    coverage_df = compute_genome_pathway_coverage(nonseed_df, pathway_to_cpds)
    print(coverage_df.head())
    print(coverage_df.tail())

    # Save for downstream/visualization
    coverage_df.to_csv("output/genome_pathway_coverage.tsv", sep="\t", index=False)

    # Genomes with most fully covered pathways
    full = coverage_df[coverage_df["coverage"] == 1]
    top_generalists = full.groupby("patric_id")["KEGG_pathway"].count().sort_values(ascending=False)
    print("\nTop 10 genomes with most fully-covered pathways:")
    print(top_generalists.head(10))

    # Top pathways covered fully by most genomes ---------
    top_pathways = full.groupby("KEGG_pathway")["patric_id"].count().sort_values(ascending=False)
    print("\nTop 10 most universally covered pathways:")
    print(top_pathways.head(10))

    # Save summary per genome or pathway
    coverage_summary = coverage_df.groupby("patric_id", group_keys=False).apply(
        lambda x: (x["coverage"] > 0.9).sum()
    ).sort_values(ascending=False)
    coverage_summary.name = "n_high_coverage_pathways"
    coverage_summary.to_csv("output/high_coverage_pathways_per_genome.csv")
    print("\nSaved summary of high-coverage pathways per genome.")


if __name__ == "__main__":
    main()    