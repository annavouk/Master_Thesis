"""
Generate a comprehensive compound summary.

This script creates a summary table of metabolic compounds present in the Seed and Non-Seed datasets.
For each compound (identified by its ModelSEED ID), the table includes:

- KEGG compound IDs (from 'seedId_keggId_module.tsv')
- Compound name (from the 'name' column in 'compounds.tsv')
- Associated ModelSEED reaction IDs (from the 'id' column in 'reactions.tsv')
- KEGG reaction IDs (from the 'abbreviation' column in 'reactions.tsv')
- KEGG modules (from 'seedId_keggId_module.tsv')
- KEGG pathways (by matching KEGG modules from 'module_map_pairs.tsv' and KEGG reactions from 'KEGG.pathways')
- Dataset origin (Seed only, Non-Seed only, or both)

The final table offers an integrated view of each compound’s metabolic context across multiple biological knowledgebases.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import re
from utils import load_data
from config import (
    SEEDS_PICKLE,
    NON_SEEDS_PICKLE,
    COMPOUNDS_TSV,
    REACTIONS_TSV,
    PATHWAYS_TSV,
    KEGG_MODULE_MAPPING,
    METADATA_DIR,
    MODULE_MAP
)

def build_detailed_compound_summary():
    """
    Builds the initial compound summary, linking compounds with reactions, modules and pathways.
    Returns:
        pd.DataFrame: The summary dataframe.
    """
    compound_df = load_data(COMPOUNDS_TSV, filetype='tsv')
    reactions_df = load_data(REACTIONS_TSV, filetype='tsv')
    pathways_df = load_data(PATHWAYS_TSV, filetype='tsv')
    kegg_map_df = load_data(KEGG_MODULE_MAPPING, filetype='tsv')
    kegg_map_df.columns = ['seed_id', 'kegg_id', 'module']
    module_map_df = load_data(MODULE_MAP, filetype="tsv")
    module_map_df.columns = ['module', 'pathway']
    module_map_df['module'] = module_map_df['module'].str.replace('md:', '', regex=False)

    seed_df = load_data(SEEDS_PICKLE, filetype="pickle")
    nonseed_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")
    seed_set = set(seed_df.columns)
    nonseed_set = set(nonseed_df.columns)
    all_cpds = list(seed_set | nonseed_set)

    # KEGG.pathways: build rxn abbreviation -> pathway(s) mapping
    rxn_to_pathways = {}
    for _, row in pathways_df.iterrows():
        pathway_id = str(row['Source ID'])
        rxns = str(row['Reactions']).split('|')
        for rxn in rxns:
            rxn = rxn.strip()
            if rxn:
                rxn_to_pathways.setdefault(rxn, set()).add(pathway_id)

    rows = []

    for cpd in all_cpds:
        # Name
        name = compound_df.loc[compound_df['id'] == cpd, 'name'].values
        name = name[0] if len(name) else ""

        # KEGG ID & modules
        kegg_ids = kegg_map_df[kegg_map_df['seed_id'] == cpd]['kegg_id'].unique().tolist()
        modules = kegg_map_df[kegg_map_df['seed_id'] == cpd]['module'].unique().tolist()

        # Reactions
        rel_rxns = reactions_df[reactions_df['compound_ids'].str.contains(cpd, na=False)]
        rxn_ids = rel_rxns['id'].dropna().unique().tolist()
        kegg_rxns = rel_rxns['abbreviation'].dropna().unique().tolist()

        # Pathways via modules
        pathway_ids = set()
        for m in modules:
            pws = module_map_df[module_map_df['module'] == m]['pathway'].dropna().unique().tolist()
            pathway_ids.update(pws)

        # Pathways via KEGG reaction abbreviations
        for kegg_rxn in kegg_rxns:
            if kegg_rxn in rxn_to_pathways:
                pathway_ids.update(rxn_to_pathways[kegg_rxn])

        pathway_ids = [p for p in set(pathway_ids) if p and str(p).strip()]

        # Dataset origin
        if cpd in seed_set and cpd in nonseed_set:
            dataset = "seed/non-seed"
        elif cpd in seed_set:
            dataset = "seed"
        else:
            dataset = "non-seed"

        rows.append({
            "SEED_ID": cpd,
            "KEGG_ID": ",".join(kegg_ids),
            "compound_name": name,
            "reactions": ",".join(rxn_ids),
            "KEGG_reactions": ",".join(kegg_rxns),
            "KEGG_modules": ",".join(modules),
            "KEGG_pathways": ",".join(pathway_ids),
            "dataset": dataset,
        })

    out_path = METADATA_DIR / "compound_summary.tsv"
    cpd_df = pd.DataFrame(rows)
    cpd_df.to_csv(out_path, sep="\t", index=False)
    print(f"Saved: {out_path}")
    return cpd_df


def report_compounds_without_pathways(cpd_df):
    """
    Prints compounds that do not have any KEGG pathways assigned.
    """
    no_pathway_df = cpd_df[cpd_df["KEGG_pathways"].isna() | (cpd_df["KEGG_pathways"].str.strip() == "")]
    print(f"Compounds without KEGG pathways: {len(no_pathway_df)}")
    print("SEED_IDs:")
    print("\n".join(no_pathway_df["SEED_ID"].tolist()))


if __name__ == "__main__":
    cpd_df = build_detailed_compound_summary()
    report_compounds_without_pathways(cpd_df)

