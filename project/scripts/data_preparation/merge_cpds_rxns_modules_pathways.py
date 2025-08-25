"""
Build compound summary by integrating ModelSEED and KEGG data.

- Labels compounds as seed, non-seed, or both
- Extracts KEGG IDs from ModelSEED aliases
- Joins KEGG modules, reactions, and pathways from KEGG JSON
- Adds BRITE ontology categories (selected roots only)
- Matches KEGG reactions of interest to ModelSEED reactions
- Outputs multiple summary tables:
    * cpds_KEGG.tsv — basic mapping of ModelSEED to KEGG IDs for compounds
    * cpds_KEGG_dataset.tsv — adds dataset origin
    * compounds_summary.tsv — enriched with KEGG info + ontology
    * unclassified.tsv — compounds without ontology classification
    * reactions_of_interest.tsv — ModelSEED reactions linked to KEGG reactions
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import re

from config import (
    SEEDS_PICKLE,
    NON_SEEDS_PICKLE,
    COMPOUNDS_TSV,
    KEGG_DATA,
    UNIQUE_REACTIONS,
    REACTIONS_TSV,
    METADATA_DIR,
)
from utils import load_data
from compounds_overview import get_compound_sets


# ======================
# Helpers
# ======================
def assign_dataset_origin(cpds, seed_set, nonseed_set):
    """Label each compound based on whether it appears in seed, non-seed or both."""
    rows = []
    for cpd in cpds:
        if cpd in seed_set and cpd in nonseed_set:
            dataset = "seed/non-seed"
        elif cpd in seed_set:
            dataset = "seed"
        elif cpd in nonseed_set:
            dataset = "non-seed"
        else:
            dataset = "unknown"
        rows.append({"SEED_ID": cpd, "dataset": dataset})
    return pd.DataFrame(rows)


def extract_kegg_ids(compound_df, column="aliases"):
    """Matches KEGG IDs from a specified column (e.g 'aliases' from compounds.tsv) to ModelSEED IDs."""
    return (
        compound_df[column]
        .fillna("")
        .apply(
            lambda x: (
                re.findall(r"KEGG:\s*([A-Za-z0-9]+)", x)[0]
                if re.findall(r"KEGG:\s*([A-Za-z0-9]+)", x)
                else None
            )
        )
    )


def extract_all_kegg_reaction_ids(compound_df, column="aliases"):
    """Extract all KEGG reaction IDs from the given column and join them with ','. """
    return (
        compound_df[column]
        .fillna("")
        .apply(lambda x: ",".join(re.findall(r"R\d{5}", x)) if "KEGG:" in x else None)
    )


def extract_kegg_info(kegg_data):
    """Extracts KEGG modules, reactions, and pathways from KEGG JSON data."""
    rows = []
    for kegg_id, entry in kegg_data.items():
        modules = [m.split()[0] for m in entry.get("MODULE", [])]
        raw_reactions = entry.get("REACTION", [])
        reactions = []
        for r in raw_reactions:
            reactions.extend(r.strip().split())
        pathways = [p.split()[0] for p in entry.get("PATHWAY", [])]
        names = entry.get("NAME", [])
        formula = entry.get("FORMULA", [""])[0]

        rows.append(
            {
                "KEGG_ID": kegg_id,
                "KEGG_modules": ", ".join(modules),
                "KEGG_reactions": ", ".join(reactions),
                "KEGG_pathways": ", ".join(pathways),
            }
        )
    return pd.DataFrame(rows)


def extract_brite_ontology(kegg_data):
    """
    Extract BRITE ontology categories for each KEGG compound.

    Returns a dict mapping KEGG_ID -> {"ontology_primary": str, "ontology_all": str}
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

    result = {}

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
            result[kegg_id] = {
                "ontology_primary": ontologies[0],
                "ontology_all": "; ".join(ontologies),
            }
        else:
            result[kegg_id] = {
                "ontology_primary": "Unclassified",
                "ontology_all": "Unclassified",
            }

    return result


def match_kegg_reactions_to_modelseed(unique_reactions_df, modelseed_reactions_df):
    """Matches KEGG reactions of interest to ModelSEED reactions based on KEGG reaction IDs."""
    kegg_reaction_ids = set(
        unique_reactions_df["KEGG_reaction_ID"].dropna().astype(str).str.strip()
    )

    matched = modelseed_reactions_df[
        modelseed_reactions_df["KEGG_reaction_ID"].isin(kegg_reaction_ids)
    ]

    selected_columns = [
        "id",
        "abbreviation",
        "name",
        "equation",
        "reversibility",
        "direction",
        "compound_ids",
        "aliases",
    ]
    return matched[selected_columns]


def build_cpds_kegg(compound_df, cpds):
    """Return mapping SEED_ID -> KEGG_ID for given compounds."""
    mask = compound_df["id"].isin(cpds)
    compound_df.loc[mask, "KEGG_ID"] = extract_kegg_ids(
        compound_df.loc[mask], column="aliases"
    )
    return compound_df.loc[mask, ["id", "name", "KEGG_ID"]].rename(
        columns={"id": "SEED_ID", "name": "compound_name"}
    )


def build_dataset_origin(cpds, seed_set, nonseed_set):
    """Return dataset origin labels for given compounds."""
    return assign_dataset_origin(cpds, seed_set, nonseed_set)


def build_kegg_summary(summary_df, kegg_df, ontology_map):
    """Enrich summary DataFrame with KEGG info and ontology."""
    final_df = summary_df.merge(kegg_df, on="KEGG_ID", how="left")
    final_df["ontology_primary"] = final_df["KEGG_ID"].map(
        lambda x: ontology_map.get(x, {}).get("ontology_primary", "Unclassified")
    )
    final_df["ontology_all"] = final_df["KEGG_ID"].map(
        lambda x: ontology_map.get(x, {}).get("ontology_all", "Unclassified")
    )
    return final_df


def save_tables(new_df, dataset_df, final_df, matched_reactions):
    """Save all output tables to METADATA_DIR."""
    # Basic mapping SEED -> KEGG
    new_df.to_csv(
        METADATA_DIR / "cpds_KEGG.tsv",
        sep="\t",
        index=False,
    )

    # Add dataset origin
    pd.merge(new_df, dataset_df, on="SEED_ID", how="left").to_csv(
        METADATA_DIR / "cpds_KEGG_dataset.tsv",
        sep="\t",
        index=False,
    )

    # Full enriched summary
    final_df.to_csv(
        METADATA_DIR / "compounds_summary.tsv",
        sep="\t",
        index=False,
    )

    # Only unclassified compounds
    final_df[final_df["ontology_primary"] == "Unclassified"].to_csv(
        METADATA_DIR / "unclassified.tsv",
        sep="\t",
        index=False,
    )

    # Reactions of interest
    matched_reactions.to_csv(
        METADATA_DIR / "reactions_of_interest.tsv",
        sep="\t",
        index=False,
    )


# ======================
# Load Inputs
# ======================
def load_inputs():
    """Load all required input files."""
    seed_df = load_data(SEEDS_PICKLE, filetype="pickle")
    nonseed_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")
    compound_df = load_data(COMPOUNDS_TSV, filetype="tsv")
    reactions_df = load_data(REACTIONS_TSV, filetype="tsv")
    kegg_data = load_data(KEGG_DATA, filetype="json")
    unique_reactions = load_data(UNIQUE_REACTIONS, filetype="tsv")
    return seed_df, nonseed_df, compound_df, reactions_df, kegg_data, unique_reactions


# ======================
# Build Compound Summary
# ======================
def prepare_summary(seed_df, nonseed_df, compound_df, kegg_data):
    """Prepare compound summary enriched with KEGG and ontology info."""
    seed_set, nonseed_set = get_compound_sets(seed_df, nonseed_df)
    cpds = list(seed_set | nonseed_set)

    new_df = build_cpds_kegg(compound_df, cpds)
    dataset_df = build_dataset_origin(cpds, seed_set, nonseed_set)
    kegg_df = extract_kegg_info(kegg_data)
    ontology_map = extract_brite_ontology(kegg_data)

    summary_df = pd.merge(new_df, dataset_df, on="SEED_ID", how="left")
    final_df = build_kegg_summary(summary_df, kegg_df, ontology_map)
    return new_df, dataset_df, final_df


# ======================
# Reaction Matching
# ======================
def prepare_reactions(reactions_df, unique_reactions):
    """Match KEGG reactions of interest to ModelSEED reactions."""
    # Extract all KEGG reaction IDs, joined with commas
    reactions_df["KEGG_reaction_IDs"] = extract_all_kegg_reaction_ids(
        reactions_df, column="aliases"
    )

    # Explode into one row per KEGG reaction ID
    exploded = reactions_df.assign(
        KEGG_reaction_IDs=reactions_df["KEGG_reaction_IDs"].str.split(",")
    ).explode("KEGG_reaction_IDs")

    # Clean spaces
    exploded["KEGG_reaction_IDs"] = exploded["KEGG_reaction_IDs"].str.strip()

    # Keep only matches with unique reactions list
    unique_ids = unique_reactions["KEGG_reaction_ID"].astype(str).str.strip()
    matched = exploded[exploded["KEGG_reaction_IDs"].isin(unique_ids)]

    # Final subset
    selected_columns = [
        "id", "abbreviation", "name", "equation", 
        "reversibility", "direction", "compound_ids", 
        "aliases", "KEGG_reaction_IDs"
    ]
    return matched[selected_columns]


# ======================
# Main
# ======================
def main():
    # Load
    seed_df, nonseed_df, compound_df, reactions_df, kegg_data, unique_reactions = (
        load_inputs()
    )

    # Build compound summary
    new_df, dataset_df, final_df = prepare_summary(
        seed_df, nonseed_df, compound_df, kegg_data
    )

    # Match reactions
    matched_reactions = prepare_reactions(reactions_df, unique_reactions)

    # Save all outputs
    save_tables(new_df, dataset_df, final_df, matched_reactions)


if __name__ == "__main__":
    main()
