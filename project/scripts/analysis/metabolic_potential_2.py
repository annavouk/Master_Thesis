"""
Metabolic Potential - Genome-by-Pathway Functional Coverage Analysis (Approach 2)

For each genome, calculates the coverage (% of compounds present as non-seeds) for each KEGG pathway.
"""


import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
from collections import defaultdict
import re

from config import (
    SEEDS_SET, 
    NON_SEEDS_SET, 
    COMPACT_METADATA, 
    COMPOUND_SUMMARY_TSV, 
    CPD_PATHWAY_MAPPING, 
    CPD_MODULE_MAPPING, 
    REACTIONS_OF_INTEREST, 
    OUTPUT_DIR,
)
from utils import load_data


# ------------------------
# Identify substrates/products of reactions of interest
# ------------------------
def parse_reaction_equation(equation, direction):
    """Parse a ModelSEED reaction equation and return substrates and products."""
    if pd.isna(equation):
        return set(), set()

    # Split into sides
    if "<=>" in equation:
        left, right = equation.split("<=>")
    elif "=>" in equation:
        left, right = equation.split("=>")
    elif "<=" in equation:  # just in case
        left, right = equation.split("<=")
    else:
        return set(), set()

    # Extract compounds
    extract = lambda side: set(re.findall(r"(cpd\d+)", side))
    left_cpds, right_cpds = extract(left), extract(right)

    # Use direction to orient
    if direction == "=":   # reversible
        substrates, products = left_cpds | right_cpds, left_cpds | right_cpds
    elif direction == ">": # forward
        substrates, products = left_cpds, right_cpds
    elif direction == "<": # backward
        substrates, products = right_cpds, left_cpds
    else:                  # unknown
        substrates, products = left_cpds, right_cpds

    return substrates, products


# ------------------------
# Genome level reaction coverage
# ------------------------
def compute_reaction_coverage(seeds_df, nonseeds_df, reactions_df):
    """Compute reaction coverage per genome."""
    results = []

    for genome in seeds_df.index:
        seed_set = set(seeds_df.loc[genome, "SeedSet"])
        nonseed_set = set(nonseeds_df.loc[genome, "NonSeedSet"])

        for _, row in reactions_df.iterrows():
            substrates = row["substrates"]

            if isinstance(substrates, str):
                substrates = set(substrates.split(";")) if substrates else set()
            elif not isinstance(substrates, set):
                substrates = set()

            if not substrates:
                continue

            n_total = len(substrates)
            n_nonseed = len(substrates & nonseed_set)
            n_seed = len(substrates & seed_set)
            coverage = n_nonseed / n_total if n_total else 0

            results.append((genome, row["KEGG_reaction_IDs"], coverage, n_total, n_nonseed, n_seed))

    return pd.DataFrame(results, columns=["genome", "reaction", "coverage", "n_total", "n_nonseed", "n_seed"])


#
def build_seed2kegg(summary_df):
    """Return dict SEED_ID -> KEGG_ID from compounds_summary_1.tsv dataframe."""
    df = summary_df.dropna(subset=["SEED_ID","KEGG_ID"]).copy()
    df["SEED_ID"] = df["SEED_ID"].astype(str).str.strip()
    df["KEGG_ID"] = df["KEGG_ID"].astype(str).str.strip()
    return (df.drop_duplicates(subset=["SEED_ID"])
              .set_index("SEED_ID")["KEGG_ID"]
              .to_dict())


def build_group_to_kegg(mapping_df: pd.DataFrame, group_col="group", compound_col="compound") :
    """
    Return dict group_id -> set(KEGG_ID) από df με δύο στήλες:
      group_col: π.χ. 'KEGG_pathway' ή 'KEGG_module' (μπορεί να περιέχει prefixes όπως 'path:')
      compound_col: π.χ. 'compound' τύπου 'cpd:C00022'
    """
    df = mapping_df.copy()
    df[group_col] = df[group_col].astype(str).str.split(":", n=1).str[-1].str.strip()
    df["KEGG_ID"] = df[compound_col].astype(str).str.replace("cpd:", "", regex=False).str.strip()
    df = df[(df[group_col] != "") & (df["KEGG_ID"] != "")]
    return df.groupby(group_col)["KEGG_ID"].agg(set).to_dict()


def compute_currency_from_groups(group_to_kegg: dict[str, set], frac_threshold: float,
                                 min_abs_groups: int | None = None) -> set:
    """Currency = compounds appearing in ≥ frac_threshold of groups (and ≥ min_abs_groups if given)."""
    if not group_to_kegg:
        return set()
    counts = defaultdict(int)
    for s in group_to_kegg.values():
        for c in s:
            counts[c] += 1
    n_groups = len(group_to_kegg)
    out = set()
    for c, k in counts.items():
        if (k / n_groups) >= frac_threshold and (min_abs_groups is None or k >= min_abs_groups):
            out.add(c)
    return out


def subtract_currency(group_to_kegg: dict[str, set], currency: set) -> dict[str, set]:
    """Επιστρέφει νέο dict χωρίς τα currency compounds."""
    if not currency:
        return dict(group_to_kegg)
    return {g: (s - currency) for g, s in group_to_kegg.items() if len(s - currency) > 0}


def genome_kegg_sets_from_matrices(seeds_df: pd.DataFrame, non_seeds_df: pd.DataFrame,
                                   seed2kegg: dict) -> tuple[dict, dict]:
    """Map 1s (SEED_ID) ανά genome σε KEGG_ID sets (returns: genome_non_kegg, genome_seed_kegg)."""
    def row_to_kegg_set(row):
        seed_ids = set(row.index[row.values == 1])
        return {seed2kegg[s] for s in seed_ids if s in seed2kegg}
    genome_non = {g: row_to_kegg_set(non_seeds_df.loc[g]) for g in non_seeds_df.index}
    genome_seed = {g: row_to_kegg_set(seeds_df.loc[g])      for g in seeds_df.index}
    return genome_non, genome_seed


def coverage_counts(genome_non_kegg: dict, genome_seed_kegg: dict,
                    group_to_kegg: dict[str, set], group_label: str,
                    min_pct_nonseed: float | None = None) -> pd.DataFrame:
    """
    Επιστρέφει long table:
      patric_id | group_label | total_in_group | n_nonseed | n_seed | pct_nonseed | pct_seed
    Δεν γίνεται κανένα I/O εδώ.
    """
    rows = []
    for genome, non_k in genome_non_kegg.items():
        seed_k = genome_seed_kegg.get(genome, set())
        for gid, grp in group_to_kegg.items():
            total = len(grp)
            if total == 0:
                continue
            n_non  = len(non_k & grp)
            n_seed = len(seed_k & grp)
            pct_non = n_non / total
            if (min_pct_nonseed is not None) and (pct_non < min_pct_nonseed):
                continue
            rows.append({
                "patric_id": genome,
                group_label: gid,
                "total_in_group": total,
                "n_nonseed": n_non,
                "n_seed": n_seed,
                "pct_nonseed": pct_non,
                "pct_seed": n_seed / total
            })
    return pd.DataFrame(rows)


def pivot_matrix(df_long: pd.DataFrame, index_col: str, col_col: str, value_col: str) -> pd.DataFrame:
    """Pivot helper (no I/O)."""
    return df_long.pivot(index=index_col, columns=col_col, values=value_col)


def compute_genome_pathway_coverage(nonseed_df, pathway_to_kegg_cpds, kegg_to_seed):
    results = []
    for genome in nonseed_df.index:
        genome_cpds = set(nonseed_df.columns[nonseed_df.loc[genome] == 1])
        for pw, kegg_cpds in pathway_to_kegg_cpds.items():
            seed_ids_in_pw = {kegg_to_seed.get(kegg_id) for kegg_id in kegg_cpds if kegg_to_seed.get(kegg_id)}
            n_required = len(kegg_cpds)  # βάση KEGG compounds
            n_covered = len(genome_cpds & seed_ids_in_pw)
            frac = n_covered / n_required if n_required else 0
            results.append({
                "patric_id": genome,
                "KEGG_pathway": pw,
                "n_required": n_required,
                "n_covered": n_covered,
                "coverage": frac
            })
    return pd.DataFrame(results)


# ------------------------
# Main
# ------------------------
def main():
    outdir = Path(OUTPUT_DIR); outdir.mkdir(parents=True, exist_ok=True)

    # Load data 
    seeds_df = load_data(SEEDS_SET, filetype="pickle")
    non_seeds_df = load_data(NON_SEEDS_SET, filetype="pickle")
    summary_df = pd.read_csv(COMPOUND_SUMMARY_TSV, sep="\t")
    pw_map_df = pd.read_csv(CPD_PATHWAY_MAPPING, sep="\t", header=None, names=["KEGG_pathway","compound"])
    mod_map_df = pd.read_csv(CPD_MODULE_MAPPING,  sep="\t", header=None, names=["KEGG_module","compound"])
    roi_df = load_data(REACTIONS_OF_INTEREST, filetype="tsv")
    meta_df = pd.read_csv(COMPACT_METADATA, dtype={"patric_id": str}, low_memory=False)

    # Identify substrates/products compounds of reactions of interest
    roi_df[["substrates", "products"]] = roi_df.apply(
    lambda row: pd.Series(parse_reaction_equation(row["equation"], row["direction"])),
    axis=1
)
    roi_df["substrates_str"] = roi_df["substrates"].apply(lambda s: ";".join(sorted(s)))
    roi_df["products_str"]   = roi_df["products"].apply(lambda s: ";".join(sorted(s)))
    print(roi_df.head()[["KEGG_reaction_IDs","substrates_str","products_str"]])

    # Reaction coverage per genome 
    df_rxn = compute_reaction_coverage(
    seeds_df.head(3),
    non_seeds_df.head(3),
    roi_df.head(20)
)
    print(df_rxn.head())


    # Build mappings in-memory
    seed2kegg  = build_seed2kegg(summary_df)
    print(f"Mapped {len(seed2kegg)} SEED_IDs to KEGG_IDs.")
    pw2kegg    = build_group_to_kegg(pw_map_df, group_col="KEGG_pathway", compound_col="compound")
    print(f"Loaded {len(pw2kegg)} KEGG pathways.")
    mod2kegg   = build_group_to_kegg(mod_map_df, group_col="KEGG_module",  compound_col="compound")
    print(f"Loaded {len(mod2kegg)} KEGG modules.")

    # Data-driven currency (ρυθμίζεις thresholds εδώ)
    currency   = compute_currency_from_groups(pw2kegg, frac_threshold=0.5, min_abs_groups=10) \
               | compute_currency_from_groups(mod2kegg, frac_threshold=0.5, min_abs_groups=5)

    pw2kegg_f  = subtract_currency(pw2kegg, currency)
    mod2kegg_f = subtract_currency(mod2kegg, currency)

    # --- Genome KEGG sets ---
    genome_non_kegg, genome_seed_kegg = genome_kegg_sets_from_matrices(seeds_df, non_seeds_df, seed2kegg)

    # --- Coverage (only data-in/data-out) ---
    df_pw  = coverage_counts(genome_non_kegg, genome_seed_kegg, pw2kegg_f,  "KEGG_pathway", min_pct_nonseed=None)
    df_mod = coverage_counts(genome_non_kegg, genome_seed_kegg, mod2kegg_f, "KEGG_module",  min_pct_nonseed=None)

    # Optional merge taxonomy (εκτός helpers)
    if meta_df is not None and not meta_df.empty:
        meta_small = meta_df[["patric_id","gtdb_taxonomy","genome_length"]].copy() if "gtdb_taxonomy" in meta_df.columns else meta_df
        df_pw  = df_pw.merge(meta_small,  on="patric_id", how="left")
        df_mod = df_mod.merge(meta_small, on="patric_id", how="left")

    # --- Save (I/O έξω από helpers) ---
    df_pw.to_csv(outdir/"coverage_counts_by_pathway.tsv", sep="\t", index=False)
    df_mod.to_csv(outdir/"coverage_counts_by_module.tsv",  sep="\t", index=False)

    pivot_pw  = pivot_matrix(df_pw,  "patric_id","KEGG_pathway","pct_nonseed")
    pivot_mod = pivot_matrix(df_mod, "patric_id","KEGG_module","pct_nonseed")
    pivot_pw.to_csv(outdir/"coverage_pct_nonseed_by_pathway_matrix.tsv", sep="\t")
    pivot_mod.to_csv(outdir/"coverage_pct_nonseed_by_module_matrix.tsv",  sep="\t")

    # Also save the detected currency list for transparency
    pd.Series(sorted(currency)).to_csv(outdir/"currency_compounds_auto.tsv", index=False, header=False)

if __name__ == "__main__":
    main()
