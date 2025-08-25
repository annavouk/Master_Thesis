"""
Prune GTDB tree to include only selected genomes using Biopython.
- Reads metadata CSV to get accession IDs
- Optionally filters a subset (e.g., top providers)
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
from Bio import Phylo
from config import COMPACT_METADATA, GTDB_TREES_DIR
from utils import load_data, parse_taxonomy

# ---------------- Parameters ----------------
TREE_TYPE = "bac"      # "bac" or "ar"
N_TOP = None            # or None to keep all
SUBSET_IDS_FILE = None # or path to file with genome ids to keep (one per line)
# --------------------------------------------

# Paths
METADATA_CSV = load_data(COMPACT_METADATA, filetype="csv")
BAC_TREE_PATH = GTDB_TREES_DIR / "bac120_r207.tree"
AR_TREE_PATH  = GTDB_TREES_DIR / "ar53_r207.tree"
OUTPUT_TREE_BAC = GTDB_TREES_DIR / "bac120_subset_pruned.tree"
OUTPUT_TREE_AR  = GTDB_TREES_DIR / "ar53_subset_pruned.tree"

def get_target_accessions(df, tree_type, subset_ids_file=None, n_top=None):
    """
    Returns set of genome accession IDs (RS_GCF...) to keep in the tree.
    """
    df = parse_taxonomy(df)
    target_domain = "Bacteria" if tree_type == "bac" else "Archaea"
    df = df[df["domain"] == target_domain]

    ids = df['accession'].astype(str).tolist()

    if subset_ids_file:
        with open(subset_ids_file) as f:
            subset = set(x.strip() for x in f if x.strip())
        ids = [i for i in ids if i in subset]

    if n_top:
        ids = ids[:n_top]

    return set(ids)

def prune_tree_biopython(tree, keep_ids, output_tree):
    """
    Prunes a Biopython Tree to keep only the leaves in keep_ids.
    Saves to output_tree (path).
    """
    all_leaves = {term.name for term in tree.get_terminals()}
    to_remove = [leaf for leaf in all_leaves if leaf not in keep_ids]
    for leaf in to_remove:
        try:
            tree.prune(target=leaf)
        except Exception:
            continue  # Ignore if not found or already pruned
    Phylo.write(tree, str(output_tree), "newick")
    print(f"Pruned tree written to: {output_tree}")

def main():
    if TREE_TYPE == "bac":
        tree_path = BAC_TREE_PATH
        output_tree = OUTPUT_TREE_BAC
    else:
        tree_path = AR_TREE_PATH
        output_tree = OUTPUT_TREE_AR

    # Read tree with Biopython
    print(f"Reading tree from: {tree_path}")
    tree = Phylo.read(str(tree_path), "newick")
    print(f"Loaded tree with {len(list(tree.get_terminals()))} leaves.")

    # Get genome IDs to keep
    keep_ids = get_target_accessions(METADATA_CSV, TREE_TYPE, subset_ids_file=SUBSET_IDS_FILE, n_top=N_TOP)
    print(f"Will keep {len(keep_ids)} genome IDs (matching the tree leaf names).")

    # Prune and save
    prune_tree_biopython(tree, keep_ids, output_tree)

if __name__ == "__main__":
    main()
