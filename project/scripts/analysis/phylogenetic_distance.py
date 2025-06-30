from Bio import Phylo
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skbio.stats.distance import mantel
from scipy.spatial.distance import pdist, squareform

# === 1. Load Trees ===
ar_tree_path = "/home/annavouk/master_thesis/gtdb_trees/ar53_r207.tree"
bac_tree_path = "/home/annavouk/master_thesis/gtdb_trees/bac120_r207.tree"
ar_tree = Phylo.read(ar_tree_path, "newick")
bac_tree = Phylo.read(bac_tree_path, "newick")

# === 2. Load Metadata & Metabolic ===
df = pd.read_csv("/home/annavouk/master_thesis/gtdb_trees/genome_size_merged_metadata.csv", low_memory=False)
total_seeds = pd.read_csv("/home/annavouk/master_thesis/project/output/metabolic_potential_summary_no_metadata.csv", low_memory=False)
df = df.merge(total_seeds, on='patric_id', how='left')
df['assembly_accession'] = df['assembly_accession'].astype(str)

# === 3. Detect Tree Prefixes (e.g. GB, RS) ===
def find_tree_prefixes(tree):
    all_prefixes = set()
    for leaf in tree.get_terminals():
        if "_" in leaf.name:
            prefix = leaf.name.split("_")[0]
            all_prefixes.add(prefix)
    return list(all_prefixes)

ar_prefixes = find_tree_prefixes(ar_tree)
bac_prefixes = find_tree_prefixes(bac_tree)
print("Tree prefixes detected (Archaea):", ar_prefixes)
print("Tree prefixes detected (Bacteria):", bac_prefixes)

def get_tree_label(acc, available_prefixes):
    if pd.isna(acc): return None
    # Use the first prefix from tree
    prefix = available_prefixes[0] if available_prefixes else 'GB'
    return prefix + '_' + acc

df['tree_label_ar'] = df['assembly_accession'].apply(lambda x: get_tree_label(x, ar_prefixes))
df['tree_label_bac'] = df['assembly_accession'].apply(lambda x: get_tree_label(x, bac_prefixes))

leaf_labels_ar = [leaf.name for leaf in ar_tree.get_terminals()]
leaf_labels_bac = [leaf.name for leaf in bac_tree.get_terminals()]
df_ar = df[df['tree_label_ar'].isin(leaf_labels_ar)].set_index('tree_label_ar')
df_bac = df[df['tree_label_bac'].isin(leaf_labels_bac)].set_index('tree_label_bac')

print("Matched Archaea:", df_ar.shape[0], "Matched Bacteria:", df_bac.shape[0])

def get_phylo_distance_matrix(tree, labels):
    """Cophenetic distance matrix for the given leaf labels"""
    n = len(labels)
    matrix = np.zeros((n, n))
    for i, l1 in enumerate(labels):
        for j, l2 in enumerate(labels):
            if i <= j:
                d = tree.distance(l1, l2)
                matrix[i, j] = matrix[j, i] = d
    return matrix

def mantel_and_scatter(tree, df_sub, label, color='teal'):
    ordered_labels = [l for l in df_sub.index if l in [leaf.name for leaf in tree.get_terminals()]]
    if len(ordered_labels) < 5:
        print(f"{label}: Too few genomes matched! Skipping.")
        return
    # Metabolic
    metab_vector = df_sub.loc[ordered_labels, 'Total_Seeds']
    metab_dist = pdist(metab_vector.values.reshape(-1, 1), metric='euclidean')
    metab_dist_matrix = squareform(metab_dist)
    # Phylogenetic
    phylo_dist_matrix = get_phylo_distance_matrix(tree, ordered_labels)
    # Mantel
    r, p, n = mantel(phylo_dist_matrix, metab_dist_matrix, method='pearson', permutations=999)
    print(f"{label} Mantel test R: {r:.3f}, p-value: {p:.4g}, n={n}")
    # Scatter
    phylo_flat = phylo_dist_matrix[np.triu_indices_from(phylo_dist_matrix, k=1)]
    metab_flat = metab_dist_matrix[np.triu_indices_from(metab_dist_matrix, k=1)]
    plt.figure(figsize=(8, 6))
    plt.scatter(phylo_flat, metab_flat, s=10, alpha=0.5, color=color)
    plt.xlabel("Phylogenetic distance")
    plt.ylabel("Metabolic distance (Seeds)")
    plt.title(f"{label}: Phylogenetic vs Metabolic Distance\n(Mantel R={r:.2f}, p={p:.2g})")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

print("\n[ARCHAEA]")
mantel_and_scatter(ar_tree, df_ar, "Archaea", color='royalblue')
print("\n[BACTERIA]")
mantel_and_scatter(bac_tree, df_bac, "Bacteria", color='crimson')
