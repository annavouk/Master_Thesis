def parse_taxonomy(data):
    """Extract taxonomic ranks from GTDB taxonomy strings."""
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    for r in ranks:
        pattern = rf"{r[0]}__([^;]+)"
        data[r] = data["gtdb_taxonomy"].str.extract(pattern)
    return data


def split_and_clean_taxonomy(df, taxonomy_col):
    """
    Split gtdb_taxonomy string to taxonomic levels, with each level started with e.g. d__
    for domain and seperated by ";".
    """
    taxonomy_series = df.loc[:, taxonomy_col].astype(str)
    taxonomy_split = taxonomy_series.str.split(";", expand=True)

    levels = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    taxonomy_split.columns = levels[: taxonomy_split.shape[1]]

    taxonomy_cleaned = taxonomy_split.apply(
        lambda col: col.str.replace(r"^[a-z]__", "", regex=True)
    )
    return taxonomy_cleaned
