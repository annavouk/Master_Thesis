import pandas as pd
import re
from pathlib import Path


def get_compound_sets(seed_df, nonseed_df):
    """Extract sets of compounds from each binary matrix."""
    return set(seed_df.columns), set(nonseed_df.columns)


# Define regex for KEGG compounds
KEGG_CPD_RE = re.compile(r"KEGG:\s*(C\d{5})")

def extract_kegg_ids(compound_df, column="aliases"):
    """
    Extract KEGG compound IDs (Cxxxxx) from a specified column.
    Looks explicitly for 'KEGG: Cxxxxx' patterns.
    Returns the first compound per row, or None if none found.
    """
    return (
        compound_df[column]
        .fillna("")
        .apply(
            lambda x: (matches[0] if isinstance(x, str) and (matches := KEGG_CPD_RE.findall(x)) else None)
        )
    )
