import pandas as pd


def split_and_clean_taxonomy(df, taxonomy_col):
    """
    Split gtdb_taxonomy string to taxonomic levels.

    Args:
        df (pd.Dataframe): gtdb_taxonomy column with each level started with e.g. d__
        for domain and seperated by ";".

    Returns:
        df (pd.DataFrame): A new DataFrame where each taxonomy rank occupies its own column,
        with prefixes removed.
    """
    taxonomy_series = df.loc[:, taxonomy_col].astype(str)
    taxonomy_split = taxonomy_series.str.split(';', expand=True)

    levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxonomy_split.columns = levels[:taxonomy_split.shape[1]]

    taxonomy_cleaned = taxonomy_split.apply(lambda col: col.str.replace(r'^[a-z]__', '', regex=True))
    return taxonomy_cleaned
