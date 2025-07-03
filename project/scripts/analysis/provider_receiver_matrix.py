"""
Metabolic Interaction Analysis Matrix

For every pair of genomes, computes the number of seed compounds (needs) of genome B
that can be supplied as non-seeds (produced) by genome A.
Outputs the 'feeding matrix' (rows: donors/providers, columns: recipients).
"""
import sys
from pathlib import Path    

sys.path.append(str(Path(__file__).resolve().parents[2]))

import pandas as pd
import numpy as np
from config import SEEDS_PICKLE, NON_SEEDS_PICKLE, OUTPUT_DIR
from utils import load_data

def compute_feeding_matrix(seeds_df, non_seeds_df):
    """
    Compute feeding matrix: For each (A, B), count how many seeds of B are non-seeds in A.

    Args:
        seeds_df (pd.DataFrame): genomes x compounds, 1=seed (needed).
        non_seeds_df (pd.DataFrame): genomes x compounds, 1=non-seed (produced).

    Returns:
        pd.DataFrame: feeding_matrix[A, B] = # compounds from seeds of B that A can provide.
    """
    # Intersect only common compounds
    common_cpds = list(set(seeds_df.columns) & set(non_seeds_df.columns))
    seeds_bin = seeds_df[common_cpds].values
    nonseeds_bin = non_seeds_df[common_cpds].values

    # Matrix multiply: for each donor (row), recipient (col)
    feeding_matrix = np.dot(nonseeds_bin, seeds_bin.T)
    return pd.DataFrame(feeding_matrix, index=non_seeds_df.index, columns=seeds_df.index)

def feeding_matrix_to_edgelist(feeding_df, threshold=1):
    """
    Converts feeding matrix to edge list (for graph/network export).
    Only keeps edges where feeding score >= threshold.

    Args:
        feeding_df (pd.DataFrame): Feeding matrix.
        threshold (int): Minimum score to keep an edge.

    Returns:
        pd.DataFrame: Edge list with columns [provider, receiver, score]
    """
    edges = []
    for donor in feeding_df.index:
        for recipient in feeding_df.columns:
            score = feeding_df.at[donor, recipient]
            if score >= threshold and donor != recipient:
                edges.append({'provider': donor, 'receiver': recipient, 'score': score})
    return pd.DataFrame(edges)

def plot_feeding_heatmap(feeding_df, top_n=30):
    """
    Plot a heatmap for top N providers/receivers.

    Args:
        feeding_df (pd.DataFrame): Feeding matrix.
        top_n (int): Show top N rows/columns.
    """
    # Top providers & receivers
    top_providers = feeding_df.sum(axis=1).sort_values(ascending=False).head(top_n).index
    top_receivers = feeding_df.sum(axis=0).sort_values(ascending=False).head(top_n).index
    sub = feeding_df.loc[top_providers, top_receivers]
    import seaborn as sns
    import matplotlib.pyplot as plt
    plt.figure(figsize=(12, 10))
    sns.heatmap(sub, cmap="viridis")
    plt.title("Feeding Matrix Heatmap (Top Providers/Receivers)")
    plt.xlabel("Recipient Genome")
    plt.ylabel("Provider Genome")
    plt.tight_layout()
    plt.show()

def main():
    # Load data
    seeds_df = load_data(SEEDS_PICKLE, filetype="pickle")
    non_seeds_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")
    
    # Try in subset
    seeds_df = seeds_df.head(200)
    non_seeds_df = non_seeds_df.head(200)


    # Compute feeding matrix 
    feeding_df = compute_feeding_matrix(seeds_df, non_seeds_df)
    feeding_df.to_csv(OUTPUT_DIR / "feeding_matrix.csv")

    print("Feeding matrix calculated and saved.")

    # Plot heatmap for summary
    plot_feeding_heatmap(feeding_df, top_n=30)

    # Print most 'generous' and 'needy' genomes
    print("\nTop providers (can feed most others):")
    print(feeding_df.sum(axis=1).sort_values(ascending=False).head(10))

    print("\nTop receivers (most needs met by others):")
    print(feeding_df.sum(axis=0).sort_values(ascending=False).head(10))

if __name__ == "__main__":
    main()
