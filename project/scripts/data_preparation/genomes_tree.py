"""
Representative Taxonomic Tree for Top Abundant Phyla

- Picks top-N phyla by abundance from COMPACT_METADATA
- Samples up to M genomes per phylum deterministically
- Builds a taxonomy-based Newick (no branch lengths)
- Exports iTOL-friendly leaves: last_rank + "_" + genome_id
- (Optional) Renders the tree with ete3 if available

Notes:
- Expects GTDB-style taxonomy strings (d__...;p__...;...;s__...)
- Skips rows with missing/malformed taxonomy
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[2]))

import argparse
import re
from collections import defaultdict, Counter

import pandas as pd
from config import COMPACT_METADATA, OUTPUT_DIR
from utils import load_data


def parse_args():
    p = argparse.ArgumentParser(description="Build taxonomy-based Newick for top phyla")
    p.add_argument("--top-n-phyla", type=int, default=5, help="How many top phyla to include")
    p.add_argument("--genomes-per-phylum", type=int, default=50, help="Max genomes per phylum")
    p.add_argument("--tax-col", default="gtdb_taxonomy", help="Column with GTDB taxonomy strings")
    p.add_argument("--id-col", default="patric_id", help="Genome identifier column")
    p.add_argument("--export-newick", default=str(OUTPUT_DIR / "tree_top_phyla.nwk"),
                   help="Output .nwk path (set to 'None' to skip)")
    p.add_argument("--render-png", default=None,
                   help="If set (e.g. path.png), tries to render with ete3 (if installed)")
    p.add_argument("--random-state", type=int, default=1, help="Deterministic sampling seed")
    return p.parse_args()


# ------------------------
# Data selection
# ------------------------
def extract_phylum(tax_string: str) -> str | None:
    if not isinstance(tax_string, str):
        return None
    m = re.search(r"p__([^;]+)", tax_string)
    return m.group(1) if m else None


def select_top_phyla(df: pd.DataFrame, tax_col: str, top_n: int) -> list[str]:
    tmp = df[tax_col].dropna().apply(extract_phylum)
    top = tmp.value_counts().nlargest(top_n).index.tolist()
    return top


def sample_per_phylum(df: pd.DataFrame, tax_col: str, id_col: str,
                      top_phyla: list[str], per_phylum: int, seed: int) -> pd.DataFrame:
    df = df.copy()
    df["phylum"] = df[tax_col].apply(extract_phylum)
    df = df[df["phylum"].isin(top_phyla)].dropna(subset=[tax_col, id_col])
    # deterministic sampling per group
    out = (
        df.groupby("phylum", group_keys=False)
          .apply(lambda g: g.sample(n=min(per_phylum, len(g)), random_state=seed))
          .reset_index(drop=True)
    )
    return out


# ------------------------
# Taxonomy -> paths -> tree
# ------------------------
def taxonomy_to_path(tax_str: str, leaf_label: str) -> list[str]:
    """
    Returns list of ranks ending with the provided leaf label.
    Filters empty tokens, trims trailing ';'
    """
    tax = (tax_str or "").rstrip(";")
    ranks = [t.strip() for t in tax.split(";") if t.strip()]
    if not ranks:
        return []
    return ranks[:-1] + [leaf_label]  # replace species token with leaf label


def build_prefix_tree(paths: list[list[str]]) -> dict:
    root: dict = {}
    for path in paths:
        d = root
        for node in path:
            d = d.setdefault(node, {})
    return root


def dict_to_newick(d: dict, name: str | None = None) -> str:
    """
    Correct Newick:
      - leaves: just the name
      - internal: '(' + comma-joined(children) + ')' + name_if_any
    """
    if not d:
        # leaf
        return (name or "")
    children = []
    for child_name, child_dict in d.items():
        children.append(dict_to_newick(child_dict, child_name))
    newick_children = ",".join(children)
    return f"({newick_children}){name or ''}"


# ------------------------
# Utilities
# ------------------------
def make_leaf_label(last_rank: str, genome_id: str) -> str:
    # iTOL-friendly, avoids spaces
    safe_id = str(genome_id).replace(" ", "_")
    return f"{last_rank}_{safe_id}"


def ensure_unique_labels(labels: list[str]) -> list[str]:
    cnt = Counter(labels)
    if all(v == 1 for v in cnt.values()):
        return labels
    seen = defaultdict(int)
    unique = []
    for lab in labels:
        seen[lab] += 1
        unique.append(lab if cnt[lab] == 1 else f"{lab}__dup{seen[lab]}")
    return unique


# ------------------------
# Main
# ------------------------
def main():
    args = parse_args()

    df = load_data(COMPACT_METADATA, filetype="csv")

    # Select & sample
    top_phyla = select_top_phyla(df, args.tax_col, args.top_n_phyla)
    print("Selected top phyla:", top_phyla)
    df_sub = sample_per_phylum(
        df, args.tax_col, args.id_col, top_phyla, args.genomes_per_phylum, args.random_state
    )
    print(f"Total genomes in subset: {len(df_sub)}")

    # Build paths
    # 1) compute provisional leaf labels (last rank + id)
    last_ranks = df_sub[args.tax_col].str.rstrip(";").str.split(";").str[-1]
    labels = [make_leaf_label(r, gid) for r, gid in zip(last_ranks, df_sub[args.id_col])]
    labels = ensure_unique_labels(labels)

    paths = []
    for tax, leaf in zip(df_sub[args.tax_col], labels):
        p = taxonomy_to_path(tax, leaf)
        if p:
            paths.append(p)

    # Tree → Newick
    tree = build_prefix_tree(paths)
    newick = dict_to_newick(tree) + ";"

    # Export
    if args.export_newick and args.export_newick.lower() != "none":
        out_path = Path(args.export_newick)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(newick)
        print(f"Newick written to: {out_path} (load directly in iTOL)")

    # Optional render with ete3 if available
    if args.render_png:
        try:
            from ete3 import Tree, TreeStyle
            t = Tree(newick, format=1)
            ts = TreeStyle()
            ts.show_leaf_name = True
            png_path = Path(args.render_png)
            png_path.parent.mkdir(parents=True, exist_ok=True)
            t.render(str(png_path), tree_style=ts, w=2000)  # wide enough for many leaves
            print(f"PNG rendered via ete3: {png_path}")
        except Exception as e:
            print(f"[WARN] ete3 render failed: {e}\n(Proceeding without PNG.)")

    # Basic sanity checks printed
    n_leaves = newick.count(",") + 1 if "(" in newick else 1
    print(f"Approx. leaves in Newick: ~{n_leaves}")


if __name__ == "__main__":
    main()
