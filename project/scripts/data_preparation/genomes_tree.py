"""
Representative Taxonomic Tree for Top Abundant Phyla

This script selects the top N most abundant phyla from the genome metadata,
randomly samples up to M genomes per phylum, constructs a taxonomy-based
(Newick format) tree using GTDB taxonomy strings, and visualizes it with ete3.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
from config import COMPACT_METADATA, OUTPUT_DIR
from utils import load_data

# PARAMETERS
TOP_N_PHYLA = 5         # Top phyla to include
GENOMES_PER_PHYLA = 50  # Max genomes per phylum
EXPORT_NEWICK = OUTPUT_DIR / "tree_top_phyla.nwk"  # Set to None to skip export


def taxonomy_to_newick_for_itol(df, tax_col="gtdb_taxonomy", id_col="patric_id"):
    """
    Returns a Newick string with leaves like s__species_genomeid.
    """
    taxonomy_list = []
    for _, row in df.iterrows():
        tax = row[tax_col]
        genome = row[id_col]
        tax = tax.rstrip(";")
        tax_ranks = [x.strip() for x in tax.split(";") if x.strip()]
        # last rank + id
        leaf_label = f"{tax_ranks[-1]}_{genome}"
        taxonomy_list.append(tax_ranks[:-1] + [leaf_label])
    # Build tree
    def build_tree(paths):
        root = {}
        for path in paths:
            d = root
            for node in path:
                d = d.setdefault(node, {})
        return root
    def dict_to_newick(d):
        if not d:
            return ""
        if len(d) == 1 and not list(d.values())[0]:
            # leaf node
            k = list(d.keys())[0]
            return k
        return "(" + ",".join(f"{k}{dict_to_newick(v)}" for k, v in d.items()) + ")"
    paths = taxonomy_list
    tree_dict = build_tree(paths)
    newick = dict_to_newick(tree_dict) + ";"
    return newick


def main():
    df = load_data(COMPACT_METADATA, filetype='csv')
    df["phylum"] = df["gtdb_taxonomy"].str.extract(r"p__([^;]+)")

    # Find top abundant phyla
    top_phyla = df["phylum"].value_counts().nlargest(TOP_N_PHYLA).index.tolist()
    print("Selected top phyla:", top_phyla)

    # For each phylum, sample up to GENOMES_PER_PHYLA genomes
    df_top = (
        df[df["phylum"].isin(top_phyla)]
        .groupby("phylum")
        .apply(lambda x: x.sample(n=min(GENOMES_PER_PHYLA, len(x)), random_state=1), include_groups=False)
        .reset_index(drop=True)
    )

    print(f"Total genomes in subset: {df_top.shape[0]}")

    # Write newick for iTOL (species+id leaves)
    newick_str = taxonomy_to_newick_for_itol(df_top, tax_col="gtdb_taxonomy", id_col="patric_id")
    if EXPORT_NEWICK:
        with open(EXPORT_NEWICK, "w") as f:
            f.write(newick_str)
        print(f"Newick tree written to {EXPORT_NEWICK} (ready for iTOL!)")


if __name__ == "__main__":
    main()
