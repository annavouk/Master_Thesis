import pandas as pd


def parse_taxonomy(data):
    """
    Extract taxonomic ranks from GTDB taxonomy strings.
    """
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    for r in ranks:
        pattern = rf"{r[0]}__([^;]+)"
        data[r] = data["gtdb_taxonomy"].str.extract(pattern)
    return data
