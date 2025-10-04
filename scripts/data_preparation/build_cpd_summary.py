"""
Build compound summary by integrating ModelSEED and KEGG data.

- Labels compounds as Seed, non-Seed, or both (from pickles)
- Extracts KEGG IDs from ModelSEED aliases (from COMPOUNDS_TSV)
- Joins KEGG modules, reactions, and pathways from KEGG_DATA_JSON
- Adds BRITE ontology categories (selected roots only) from KEGG_DATA_JSON
- Outputs summary table:
    * compounds_summary.tsv - enriched with KEGG info + ontology
    * brite_unclassified.tsv - compounds without ontology classification
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import re

from config import (
    SEEDS_PICKLE,  # input for Seed compounds
    NON_SEEDS_PICKLE,  # input for non-Seed compounds
    COMPOUNDS_TSV,  # input data for ModelSEED compounds
    KEGG_DATA_JSON,  # input data for KEGG compounds
    COMPOUND_SUMMARY_TSV,  # output
    BRITE_UNCLASSIFIED_TSV,  # output
)
from utils import load_data, get_compound_sets, extract_kegg_ids


# ======================
# Helpers
# ======================
def assign_dataset_origin(cpds, seed_set, nonseed_set):
    """Label each compound based on whether it appears in Seed, non-Seed or both."""
    cpds_rows = []
    for cpd in cpds:
        if cpd in seed_set and cpd in nonseed_set:
            dataset = "seed/non-seed"
        elif cpd in seed_set:
            dataset = "seed"
        elif cpd in nonseed_set:
            dataset = "non-seed"
        else:
            dataset = "unknown"
        cpds_rows.append({"ModelSEED_ID": cpd, "dataset": dataset})

    cpds_dataset = pd.DataFrame(cpds_rows)
    return cpds_dataset


def extract_kegg_info(kegg_data):
    """Extracts KEGG modules, reactions and pathways from KEGG_DATA_JSON."""
    kegg_rows = []
    for kegg_id, entry in kegg_data.items():
        modules = [m.split()[0] for m in entry.get("MODULE", [])]
        raw_reactions = entry.get("REACTION", [])
        reactions = []
        for r in raw_reactions:
            reactions.extend(r.strip().split())
        pathways = [p.split()[0] for p in entry.get("PATHWAY", [])]

        kegg_rows.append(
            {
                "KEGG_ID": kegg_id,
                "KEGG_modules": ", ".join(modules),
                "KEGG_reactions": ", ".join(reactions),
                "KEGG_pathways": ", ".join(pathways),
            }
        )
    kegg_info_df = pd.DataFrame(kegg_rows)
    return kegg_info_df


def extract_brite_ontology(kegg_data):
    """
    Extract BRITE ontology categories for each KEGG compound.

    - ontology_primary: first/most specific ontology label found
    - ontology_all: all ontology labels joined by ';'
    """
    relevant_brite_roots = {
        "br08001",
        "br08002",
        "br08003",
        "br08005",
        "br08006",
        "br08007",
        "br08009",
        "br08021",
    }

    brite_rows = []

    for kegg_id, entry in kegg_data.items():
        brite = entry.get("BRITE", [])
        ontologies = []

        for i, line in enumerate(brite):
            if any(f"[BR:{root}]" in line for root in relevant_brite_roots):
                j = i + 1
                last_valid_label = None
                while j < len(brite):
                    current = brite[j].strip()
                    if re.match(r"^C\d{5}\b", current):  # compound line reached
                        if last_valid_label:
                            ontologies.append(last_valid_label)
                        break
                    elif (
                        not current
                        or current.startswith("D")
                        or current.startswith("[")
                    ):
                        pass
                    else:
                        last_valid_label = current
                    j += 1

        if ontologies:
            brite_rows.append(
                {
                    "KEGG_ID": kegg_id,
                    "ontology_primary": ontologies[0],
                    "ontology_all": "; ".join(ontologies),
                }
            )
        else:
            brite_rows.append(
                {
                    "KEGG_ID": kegg_id,
                    "ontology_primary": "Unclassified",
                    "ontology_all": "Unclassified",
                }
            )

    brite_df = pd.DataFrame(
        brite_rows, columns=["KEGG_ID", "ontology_primary", "ontology_all"]
    )
    return brite_df


def build_cpds_kegg(compound_df, cpds):
    """Return mapping ModelSEED_ID -> KEGG_ID for given compounds."""
    cpd_df = compound_df.copy()

    # Filter only rows with compounds of interst
    mask = cpd_df["id"].isin(cpds)
    subset = cpd_df.loc[mask, ["id", "name", "aliases"]].copy()

    # Extract KEGG_ID
    subset["KEGG_ID"] = extract_kegg_ids(subset, column="aliases")

    # Rename
    cpds_kegg_df = subset[["id", "name", "KEGG_ID"]].rename(
        columns={"id": "ModelSEED_ID", "name": "compound_name"}
    )
    return cpds_kegg_df


def save_tables(summary_df):
    """Save output Dataframe to METADATA_DIR."""
    # Full enriched summary
    out1 = Path(COMPOUND_SUMMARY_TSV)
    out1.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(out1, sep="\t", index=False)

    # Only unclassified compounds
    out2 = Path(BRITE_UNCLASSIFIED_TSV)
    out2.parent.mkdir(parents=True, exist_ok=True)
    summary_df[summary_df["ontology_primary"] == "Unclassified"].to_csv(
        out2, sep="\t", index=False
    )


# ======================
# Load Inputs
# ======================
def load_inputs():
    """Load all required input files."""
    seed_df = load_data(SEEDS_PICKLE, filetype="pickle")
    nonseed_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")
    compound_df = load_data(COMPOUNDS_TSV, filetype="tsv")
    kegg_data = load_data(KEGG_DATA_JSON, filetype="json")
    return seed_df, nonseed_df, compound_df, kegg_data


# ======================
# Build Compound Summary
# ======================
def prepare_summary(seed_df, nonseed_df, compound_df, kegg_data):
    """
    Build compound summary Dataframe by merging:
    - ModelSEED Seed/non-Seed compounds origin and compound name
    - KEGG data (KEGG IDs, modules, reactions, pathways)
    - BRITE ontology (primary + all)
    """
    # Extract compounds
    seed_set, nonseed_set = get_compound_sets(seed_df, nonseed_df)
    cpds = sorted(seed_set | nonseed_set)

    # ModelSEED -> KEGG IDs for compounds of interest
    cpds_kegg_df = build_cpds_kegg(
        compound_df, cpds
    )  # cols: ModelSEED_ID, compound_name, KEGG_ID

    # Dataset origin (Seed / non-Seed / both)
    dataset_df = assign_dataset_origin(
        cpds, seed_set, nonseed_set
    )  # cols: SEED_ID, dataset
    dataset_df = dataset_df.rename(columns={"SEED_ID": "ModelSEED_ID"})

    # KEGG info (modules, reactions, pathways)
    kegg_info_df = extract_kegg_info(
        kegg_data
    )  # cols: KEGG_ID, KEGG_modules, KEGG_reactions, KEGG_pathways

    # BRITE ontology
    brite_df = extract_brite_ontology(
        kegg_data
    )  # cols: KEGG_ID, ontology_primary, ontology_all

    # Merge
    summary_df = (
        cpds_kegg_df.merge(dataset_df, on="ModelSEED_ID", how="left")
        .merge(kegg_info_df, on="KEGG_ID", how="left")
        .merge(brite_df, on="KEGG_ID", how="left")
    )

    # Columns order
    cols = [
        "ModelSEED_ID",
        "compound_name",
        "KEGG_ID",
        "dataset",
        "KEGG_modules",
        "KEGG_reactions",
        "KEGG_pathways",
        "ontology_primary",
        "ontology_all",
    ]
    summary_df = summary_df.reindex(columns=cols)
    summary_df = summary_df.drop_duplicates(subset=["ModelSEED_ID"]).reset_index(
        drop=True
    )
    summary_df = summary_df.fillna("")
    return summary_df


# ======================
# Main
# ======================
def main():
    # Load
    seed_df, nonseed_df, compound_df, kegg_data = load_inputs()

    # Build compound summary
    summary_df = prepare_summary(seed_df, nonseed_df, compound_df, kegg_data)

    # Save all outputs
    save_tables(summary_df)
    print(f"Compound summary file saved to {COMPOUND_SUMMARY_TSV}")


if __name__ == "__main__":
    main()
