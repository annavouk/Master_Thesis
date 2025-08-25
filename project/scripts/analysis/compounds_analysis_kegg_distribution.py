"""
Analyze which compounds, come with higher seeds/non-seeds frequency, based on binary presence.
Calculates total genomes count per seed/non-seed, combines them and appends the KEGG modules the compound is present in.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pandas as pd
from utils import load_data
from config import SEEDS_PICKLE, NON_SEEDS_PICKLE, COMPOUND_SUMMARY_TSV, OUTPUT_DIR


def genomes_per_seed(df):
    """Count the number of genomes (value == 1) per seed compound."""
    genomes = df.sum()
    most_genomes = genomes.idxmax()
    least_genomes = genomes.idxmin()
    max_count = genomes.max()
    min_count = genomes.min()

    sum_row = pd.DataFrame(genomes).T
    sum_row.index = ['Total_Genomes_per_Seed']
    updated_df = pd.concat([df, sum_row])

    return (genomes, most_genomes, max_count, least_genomes, min_count, updated_df)


def genomes_per_non_seed(df):
    """Count the number of genomes (value == 1) per non-seed compound."""
    # Count 1s per column
    genomes = df.sum()
    most_genomes = genomes.idxmax()
    least_genomes = genomes.idxmin()
    max_count = genomes.max()
    min_count = genomes.min()

    sum_row = pd.DataFrame(genomes).T
    sum_row.index = ['Total_Genomes_per_non_Seed']
    updated_df = pd.concat([df, sum_row])

    return (genomes, most_genomes, max_count, least_genomes, min_count, updated_df)


def combine_seed_and_nonseed_counts(seed_counts, non_seed_counts):
    """Combine seed and non-seed genome counts into a single DataFrame."""
    seed_counts.name = "seed_genome_count"
    non_seed_counts.name = "non_seed_genome_count"

    combined_df = pd.concat([seed_counts, non_seed_counts], axis=1).fillna(0).astype(int)
    return combined_df


def enrich_with_kegg_info(df, mapping_df, seed_id_col='ModelSEED ID'):
    """Enrich genome count per Seed/ non-Seed DataFrame with KEGG modules and pathways."""
    # Ensure column names match
    mapping_df = mapping_df.rename(columns={"SEED_ID": "ModelSEED ID"})

    # Parse to lists
    for col in ['KEGG_modules', 'KEGG_pathways']:
        if col in mapping_df.columns:
            mapping_df[col] = mapping_df[col].fillna('')
            mapping_df[col] = mapping_df[col].apply(
                lambda x: [v.strip() for v in str(x).split(',') if v.strip()] if x else []
            )
    # Group by ModelSEED ID and merge lists (flatten unique values)
    agg = mapping_df.groupby(seed_id_col).agg({
        'KEGG_modules': lambda x: sorted(set(i for sub in x for i in sub)),
        'KEGG_pathways': lambda x: sorted(set(i for sub in x for i in sub)),
        'dataset': lambda x: sorted(set(x))
    }).reset_index()

    # Merge aggregated info into main df
    enriched = pd.merge(df, agg, on=seed_id_col, how='left')
    return enriched


def flatten_column(df, id_col, list_col):
    """Flattens a dataframe with a column of lists, returning a new dataframe
    with one row per list element."""
    # Convert string-lists to real lists if needed
    df = df.copy()
    df[list_col] = df[list_col].apply(lambda x: eval(x) if isinstance(x, str) and x.startswith('[') else x)

    # Explode
    exploded = df[[id_col, list_col]].explode(list_col)
    # Drop rows with empty/null
    exploded = exploded[exploded[list_col].notnull() & (exploded[list_col] != "")]
    return exploded


def main():
    # Load
    seeds_df = load_data(SEEDS_PICKLE, filetype="pickle")
    non_seeds_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")
    mapping_df = load_data(COMPOUND_SUMMARY_TSV, filetype="tsv")

    # Analyze seeds
    (
        seed_counts,
        most_seed,
        max_seed,
        least_seed,
        min_seed,
        updated_seeds_df
    ) = genomes_per_seed(seeds_df)

    # Analyze non-seeds
    (
        non_seed_counts,
        most_non_seed,
        max_non_seed,
        least_non_seed,
        min_non_seed,
        updated_non_seeds_df
    ) = genomes_per_non_seed(non_seeds_df)

    # Combine results
    combined_counts_df = combine_seed_and_nonseed_counts(seed_counts, non_seed_counts)
    combined_counts_df.reset_index(inplace=True)
    combined_counts_df.rename(columns={"index": "ModelSEED ID"}, inplace=True)
    #print(combined_counts_df.head())

    # Enrich
    enriched_df = enrich_with_kegg_info(combined_counts_df, mapping_df)
    #print(enriched_df.head())
    #enriched_df.to_csv(OUTPUT_DIR/'genomes_per_seed_non_seed_kegg_modules_pathways.csv', index=False)

    # Count how many KEGG pathways each seed paricipate in
    enriched_df['num_kegg_pathways'] = enriched_df['KEGG_pathways'].apply(lambda x: len(set(x)) if isinstance(x, list) else 0)

    # Count genome coverage per seed
    n_genomes = seeds_df.shape[0]
    enriched_df['percent_as_seed'] = (enriched_df['seed_genome_count'] / n_genomes * 100).round(2)
    enriched_df['percent_as_non_seed'] = (enriched_df['non_seed_genome_count'] / n_genomes * 100).round(2)
    enriched_df[['ModelSEED ID', 'percent_as_seed', 'percent_as_non_seed']].to_csv(
    OUTPUT_DIR / "genome_coverage_per_compound.csv", index=False
    )

    # Keep those compounds that are actually seeds
    seed_ids = set(seeds_df.columns)
    pathways_per_seed = enriched_df[enriched_df['ModelSEED ID'].isin(seed_ids)][['ModelSEED ID', 'num_kegg_pathways']]
    pathways_per_seed.to_csv(OUTPUT_DIR / "counts_pathways_per_seed.csv", index=False)

    # Explode KEGG pathways
    flat_pathway = flatten_column(enriched_df, "ModelSEED ID", "KEGG_pathways")
    print(flat_pathway.head())

    # Count how many seeds per pathway
    counts_pathway = flat_pathway['KEGG_pathways'].value_counts().reset_index()
    counts_pathway.columns = ['KEGG_pathway', 'num_seeds']
    counts_pathway.to_csv(OUTPUT_DIR/'counts_seeds_per_pathway.csv', index=False)

    # Explode KEGG modules
    flat_module = flatten_column(enriched_df, "ModelSEED ID", "KEGG_modules")
    #print(flat_module.head())

    # Count how many seeds per module
    counts_module = flat_module['KEGG_modules'].value_counts().reset_index()
    counts_module.columns = ['KEGG_module', 'num_seeds']
    counts_module.to_csv(OUTPUT_DIR/'counts_seeds_per_module.csv', index=False)
    
    # Construct a binary participation matrix to be used for heatmap visualization (seeds × KEGG pathways)
    matrix = {}
    for _, row in enriched_df.iterrows():
        cpd = row["ModelSEED ID"]
        pws = row["KEGG_pathways"] if isinstance(row["KEGG_pathways"], list) else []
        matrix[cpd] = {p: 1 for p in pws if p}
    pd.DataFrame.from_dict(matrix, orient="index").fillna(0).astype(int).to_csv(
        OUTPUT_DIR / "seed_pathway_participation_matrix.csv"
    )



if __name__ == "__main__":
    main()
