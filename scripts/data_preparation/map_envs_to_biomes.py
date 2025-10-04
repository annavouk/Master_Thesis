"""
Process PREGO environment annotations from genome metadata and assign broad biome categories.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import re
import pandas as pd
import numpy as np

from config import (
    COMPACT_METADATA_TSV,  # input
    COMPACT_METADATA_WITH_BIOME_TSV,  # output
)

from utils import load_data


# -------------------------------
# Rules for assigning biomes
# -------------------------------
RULES_REGEX = [
    (
        r"\bsoil\b|rhizosphere|forest soil|bulk soil|paddy field|forest\b|grassland|woodland|farm|field\b|vegetated area|agricultural feature|peatland|cropland biome",
        "Soil",
    ),
    (
        r"\bsea\b|marine|ocean|brackish|lagoon|coastal|coast|harbor|beach|estuarine biome|mangrove biome",
        "Marine",
    ),
    (
        r"\bfresh\s*water\b|\bfreshwater\b|lake|pond|wetland|ground\s*water|river|stream|spring|aquifer|watercourse",
        "Freshwater",
    ),
    (
        r"sludge|waste\b|sewage|bioreactor|landfill|compost|treatment plant|gold mine|animal manure|building|mine\b|petroleum|oil\b|urban|food processing factory|dairy|cultivated environment|village biome|anthropogenic environment",
        "Engineered",
    ),
    (
        r"glacier|permafrost|hot\s*spring|hydrothermal vent|saline\b|acid\b|tundra|extreme|desert|alkaline|high temperature|ice mass|glacial feature|fumarole",
        "Extreme",
    ),
    (r"gut|intestin|oral|skin\b|host\b|hospital", "Host"),
]

RULES_REGEX = [(re.compile(p, flags=re.I), lab) for p, lab in RULES_REGEX]

PRIORITY = ["Soil", "Marine", "Freshwater", "Engineered", "Extreme", "Host", "Other"]


# -------------------------------
# Functions
# -------------------------------
def preprocess_envs(df):
    """Explode and clean the all_envs column."""
    df["all_envs"] = df["all_envs"].fillna("").str.split(";")
    df_long = df.explode("all_envs")
    df_long["all_envs"] = df_long["all_envs"].str.strip()
    df_long["all_envs"] = df_long["all_envs"].replace("", np.nan)
    return df_long


def pick_with_priority(biomes):
    for biome in PRIORITY:
        if biome in biomes.values:
            return biome
    return "Other"


def assign_biome(env):
    """Assign a broad biome category to a raw environment string."""
    if pd.isna(env):
        return "Other"
    e = env.lower()
    for pattern, biome in RULES_REGEX:
        if re.search(pattern, e):
            return biome
    return "Other"


def summarize_envs(df_long, top_n=20):
    """Summarize environment counts and percentages."""
    counts = df_long["all_envs"].value_counts(dropna=False)
    summary = (
        counts.rename("count").reset_index().rename(columns={"index": "environment"})
    )
    summary["percent"] = 100 * summary["count"] / summary["count"].sum()
    print("\nTop environments:")
    print(summary.head(top_n))
    return summary


def add_biomes(df_long):
    """Add broad biome assignments to dataframe."""
    df_long["broad_biome"] = df_long["all_envs"].apply(assign_biome)
    return df_long


def get_main_biome(df_long):
    """Assign one main broad biome per genome using priority tie-breaker."""
    chosen = df_long.groupby("patric_id")["broad_biome"].apply(pick_with_priority)
    return chosen.reset_index().rename(columns={"broad_biome": "main_biome"})


def check_unmatched(df_long, top_n=50):
    """Show the most frequent environment terms currently assigned to Other."""
    unmatched = df_long.loc[df_long["broad_biome"] == "Other", "all_envs"]
    return unmatched.value_counts().head(top_n)


# -------------------------------
# Main
# -------------------------------
def main():
    # Load metadata using project utility
    df = load_data(COMPACT_METADATA_TSV, filetype="tsv")
    df_long = preprocess_envs(df)

    print("\nExample rows:")
    print(df_long.loc[:, ["patric_id", "all_envs"]].head(10))
    print(f"\nNumber of rows after explode: {len(df_long)}")

    # Summarize environments
    env_summary = summarize_envs(df_long, top_n=30)

    # Assign biomes
    df_long = add_biomes(df_long)

    print("\nBroad biome distribution:")
    print(df_long["broad_biome"].value_counts())

    # Assign main biome per genome
    genome_biomes = get_main_biome(df_long)

    print("\nExample genome to main biome mapping:")
    print(genome_biomes.head())
    print(f"\nNumber of genomes with assigned main biome: {len(genome_biomes)}")

    # Check unmatched environments
    print("\nTop unmatched environments:")
    print(check_unmatched(df_long, top_n=50))

    # Coverage report
    total_genomes = genome_biomes.shape[0]
    n_other = (genome_biomes["main_biome"] == "Other").sum()
    coverage = 100 * (1 - n_other / total_genomes)
    print(f"\nCoverage: {coverage:.1f}% of genomes assigned to a meaningful biome")

    # Merge biome assignments back into compact metadata
    df_with_biome = df.merge(genome_biomes, on="patric_id", how="left")
    print(f"\nMetadata with biome assignments, example rows:")
    print(df_with_biome.head())

    # Save to file
    out_path = Path(COMPACT_METADATA_WITH_BIOME_TSV)
    df_with_biome.to_csv(out_path, sep="\t", index=False)

    print(f"\nSaved merged metadata with biome assignments to {out_path}")

    print("\nGenomes per biome:")
    print(genome_biomes["main_biome"].value_counts())

    return df_long, env_summary, genome_biomes


if __name__ == "__main__":
    main()
