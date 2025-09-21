"""
Phylogeny-Trait Correlation Analysis

This script tests the relationship between phylogenetic distances and 
genome-level metabolic traits (Seeds_per_Mbp, Non_Seeds_per_Mbp, Ratio_per_Mbp).
"""

import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parents[3]))

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import statsmodels.api as sm
import random
from itertools import combinations
from skbio.stats.distance import mantel
from scipy.stats import spearmanr
from skbio import TreeNode, DistanceMatrix
from Bio import Phylo
from io import StringIO
from config import BAC_TREE_PRUNED, AR_TREE_PRUNED, METABOLIC_POTENTIAL_1, OUTPUT_DIR
from utils import load_data, parse_taxonomy

TRAITS = ["Seeds_per_Mbp", "Non_Seeds_per_Mbp", "Ratio_per_Mbp"]
DOMAIN_TREES = {
    "Bacteria": Path(BAC_TREE_PRUNED),
    "Archaea": Path(AR_TREE_PRUNED)
}

def trait_distance_matrix(vec, ids):
    v = vec.values.astype(float).reshape(-1, 1)
    d = np.sqrt((v - v.T) ** 2)
    return DistanceMatrix(d, ids)

def normalize_id(x):
    return x.strip().replace(" ", "_")

def run_domain(dom, tree_path, df, plots_dir, r_threshold=0.6):
    tree = TreeNode.read(str(tree_path))

    tip_map = {t.name.replace(" ", "_"): t.name for t in tree.tips()}
    tips_norm = set(tip_map.keys())

    df = df.copy()

    for c in TRAITS:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.replace([np.inf, -np.inf], np.nan)

    df = df[df["accession_norm"].isin(tips_norm)].dropna(subset=TRAITS)
    matched = df["accession_norm"].nunique()
    print(f"[{dom}] Matched tips after filter: {matched}")
    if matched < 20:
        return []

    all_rows = []

    for phylum, sub_df in df.groupby("phylum"):
        if sub_df["accession_norm"].nunique() < 20:
            continue

        final_norm = sorted(set(sub_df["accession_norm"]) & tips_norm)
        if len(final_norm) < 20:
            continue

        final_original = [tip_map[n] for n in final_norm if n in tip_map]
        tree_tips = {t.name for t in tree.tips()}
        final_original = [t for t in final_original if t in tree_tips]

        if len(final_original) < 20:
            print(f"[{dom} | {phylum}] Skipping: only {len(final_original)} in tree after cleanup")
            continue

        print(f"[{dom} | {phylum}] {len(final_original)} genomes in tree")

        sub_tree = tree.shear(final_original)
        phylo_dm = sub_tree.tip_tip_distances()

        tree_ids_original = list(phylo_dm.ids)
        tree_ids_norm = [normalize_id(t) for t in tree_ids_original]

        common_original = []
        order_norm = []
        for orig, norm in zip(tree_ids_original, tree_ids_norm):
            if norm in sub_df["accession_norm"].values:
                common_original.append(orig)
                order_norm.append(norm)

        phylo_dm = phylo_dm.filter(common_original)
        sub_df_ordered = sub_df.set_index("accession_norm").loc[order_norm]

        for tr in TRAITS:
            vec = sub_df_ordered[tr]
            if vec.var() == 0:
                continue
            tdm = trait_distance_matrix(vec, common_original)
            r, p, n = mantel(phylo_dm, tdm, method="spearman", permutations=9999, strict=True)

            print(f"[{dom} | {phylum}] Mantel for {tr}: r={r:.3f}, p={p:.4f}, n_pairs={n}")

            all_rows.append({
                "domain": dom,
                "phylum": phylum,
                "trait": tr,
                "r": float(r),
                "p": float(p),
                "n_pairs": int(n),
                "n_genomes": int(len(sub_df_ordered)),
            })

    return all_rows

def plot_mantel_heatmap(df, r_col='r', p_col='p', phylum_col='phylum', trait_col='trait'):
    heatmap_data = df.pivot(index=phylum_col, columns=trait_col, values=r_col)
    signif = df.pivot(index=phylum_col, columns=trait_col, values=p_col) <= 0.05

    plt.figure(figsize=(12, max(6, 0.25*len(heatmap_data))))
    plt.title("Mantel test correlation (r) per phylum and trait")

    im = plt.imshow(heatmap_data, aspect='auto', cmap="viridis")
    plt.colorbar(im, label="Mantel r")

    plt.yticks(ticks=np.arange(len(heatmap_data.index)), labels=heatmap_data.index)
    plt.xticks(ticks=np.arange(len(heatmap_data.columns)), labels=heatmap_data.columns)

    for i, phylum in enumerate(heatmap_data.index):
        for j, trait in enumerate(heatmap_data.columns):
            if signif.loc[phylum, trait]:
                plt.text(j, i, "*", ha='center', va='center', color='white', fontsize=14)

    plt.tight_layout()
    plt.show()

def compute_patristic_distances(tree, allowed_names=None):
    """Return dict of distances"""
    tips = tree.get_terminals()
    if allowed_names:
        tips = [t for t in tips if t.name in allowed_names]
    dist_dict = {}
    for a, b in combinations(tips, 2):
        dist_dict[(a.name, b.name)] = tree.distance(a, b)
    return dist_dict

def scatter_phylo_vs_traitdiff_with_trend(df, dist_dict, trait, dom, ax, sample_size=None, text_pos=(0.02, 0.80), fontsize=9):
    pairs = []
    for (g1, g2), dist in dist_dict.items():
        if g1 in df.index and g2 in df.index:
            diff = abs(df.loc[g1, trait] - df.loc[g2, trait])
            pairs.append((dist, diff))
    pairs = np.array(pairs)

    if sample_size and len(pairs) > sample_size:
        idx = np.random.choice(len(pairs), sample_size, replace=False)
        pairs = pairs[idx]

    # Hexbin
    hb = ax.hexbin(pairs[:, 0], pairs[:, 1], gridsize=50, cmap="viridis", bins="log")
    ax.set_xlabel("Phylogenetic distance")
    ax.set_ylabel(f"Absolute difference in {trait.replace('_', ' ')}")

    # Spearman
    r, p = spearmanr(pairs[:, 0], pairs[:, 1])

    ax.text(text_pos[0], text_pos[1],
            f"Spearman r={r:.2f}\np={p:.2e}",
            transform=ax.transAxes, ha="left", va="top",
            fontsize=fontsize,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, linewidth=0),
            zorder=5)

    # LOESS
    loess_frac = 0.3 if len(pairs) < 5000 else 0.1
    loess_result = sm.nonparametric.lowess(pairs[:, 1], pairs[:, 0], frac=loess_frac)
    ax.plot(loess_result[:, 0], loess_result[:, 1], color="red", lw=2, zorder=4)

    return hb


def plot_tree_with_traits_biopython(tree_file, traits_df,
                                    trait_col='Ratio_per_Mbp', label_col='phylum',
                                    output_file=None):

    safe_tree = sanitize_newick(tree_file)
    tree = Phylo.read(safe_tree, "newick")

    df = traits_df.copy()
    df['accession_norm'] = df['accession'].astype(str).str.strip().str.replace(" ", "_")
    trait_map = dict(zip(df['accession_norm'], pd.to_numeric(df[trait_col], errors='coerce')))

    # Map accession to phylum
    phylum_map = dict(zip(df['accession_norm'], df[label_col]))

    v = pd.Series(trait_map).dropna()
    if v.empty:
        print(f"[plot] No valid values for {trait_col}; skipping plot.")
        return

    vmin = max(v.min(), 1e-3)
    vmax = np.percentile(v, 99)
    norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax) if vmin > 0 else mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    cmap = mpl.colormaps.get_cmap("plasma")
    name_to_terminal = {t.name: t for t in tree.get_terminals()}

    # Colour the branches
    for name, clade in name_to_terminal.items():
        val = trait_map.get(name, None)
        if pd.notna(val):
            clade.color = mpl.colors.to_hex(cmap(norm(val)), keep_alpha=False)
        else:
            clade.color = "#999999"

    # Choose one tip per phylum
    phylum_representatives = {}
    for tip in tree.get_terminals():
        phylum = phylum_map.get(tip.name)
        if phylum and phylum not in phylum_representatives:
            phylum_representatives[phylum] = tip.name

    fig = plt.figure(figsize=(12, 16))
    ax = fig.add_subplot(1, 1, 1)

    def label_func(clade):
        if clade.is_terminal() and clade.name in phylum_representatives.values():
            return phylum_map.get(clade.name, "")
        else:
            return ""

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=label_func,
        show_confidence=False
    )

    for txt in ax.texts:
        txt.set_fontsize(8)

    for line in ax.get_lines():
        line.set_linewidth(1.0)

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(sm, ax=ax)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label(trait_col, fontsize=12)

    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"[plot] Saved: {output_file}")
    else:
        plt.show()


def cluster_trait_summary(traits_df, cluster_lists, trait_col='Seeds_per_Mbp'):
    results = []
    for i, cluster in enumerate(cluster_lists):
        subset = traits_df[traits_df['accession_norm'].isin(cluster)]
        mean_val = subset[trait_col].mean()
        results.append({"cluster": i, "n_genomes": len(subset), "mean_trait": mean_val})
    return pd.DataFrame(results)

def _max_pairwise_patristic(tree: Phylo.BaseTree.Tree, tips) -> float:
    # tips: list of Clade terminals
    if len(tips) < 2:
        return 0.0
    return max(tree.distance(a, b) for a, b in combinations(tips, 2))

def clusters_by_maxdist(tree_file: str, cutoff: float) -> list[list[str]]:
    """
    Return monophyletic clusters of tip names. A cluster is any internal clade
    whose max pairwise patristic distance among its tips <= cutoff.
    We choose the largest such clades (don’t descend further once condition holds).
    """
    safe_tree = sanitize_newick(tree_file)
    tree = Phylo.read(safe_tree, "newick")

    clusters = []

    def visit(clade):
        tips = clade.get_terminals()
        if len(tips) < 2:
            return
        maxd = _max_pairwise_patristic(tree, tips)
        if maxd <= cutoff:
            clusters.append([t.name for t in tips])
        else:
            for child in clade.clades:
                visit(child)

    visit(tree.root)
    return clusters


def main():
    outdir = Path(OUTPUT_DIR) / "plots" / "phylogeny" / "phylo_scatter"
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_data(METABOLIC_POTENTIAL_1, filetype="csv")
    df['accession_norm'] = df['accession'].astype(str).str.strip().str.replace(" ", "_")
    df = parse_taxonomy(df)
    df = df.set_index("accession_norm")

    SAMPLE_TIPS = 100

    for dom, tree_path in DOMAIN_TREES.items():
        tree = Phylo.read(str(tree_path), "newick")
        allowed = [t.name for t in tree.get_terminals() if t.name in df.index]
        if len(allowed) > SAMPLE_TIPS:
            allowed = random.sample(allowed, SAMPLE_TIPS)

        dist_dict = compute_patristic_distances(tree, allowed_names=set(allowed))

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))
        plt.suptitle(f"Phylogenetic distance vs trait difference — {dom}")

        hbs = []
        for ax, trait in zip(axes, TRAITS):
            hb = scatter_phylo_vs_traitdiff_with_trend(df, dist_dict, trait, dom, ax=ax, sample_size=50000, text_pos=(0.02, 0.78))
            hbs.append(hb)

        # Οριζόντιο colorbar κάτω από όλα
        cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.03])  # left, bottom, width, height
        cbar = fig.colorbar(hbs[0], cax=cbar_ax, orientation="horizontal")
        cbar.set_label("log10(#pairs)")

        plt.tight_layout(rect=[0, 0.08, 1, 0.95])  # αφήνει χώρο κάτω
        outpath = outdir / f"{dom.lower()}_phylo_vs_traitdiff_trend.png"
        plt.savefig(outpath, dpi=300)
        plt.close()
        print(f"[plot] Saved: {outpath}")


if __name__ == "__main__":
    raise SystemExit(main())

