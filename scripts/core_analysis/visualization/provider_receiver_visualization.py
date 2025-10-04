"""
Metabolic Interaction Visualization

This script processes edge lists of microbial interactions, computes top providers/receivers, and
plots heatmaps.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from config import (
    OUTLIERS_L2H_EDGES_CSV,  # input
    OUTLIERS_H2L_EDGES_CSV,  # input
    CYANO_HALO_EDGES_CSV,  # input
    METABOLIC_POTENTIAL_TSV,  # metadata
    INTERACTION_ANALYSIS_PLOTS_DIR,  # plots dir
)
from utils import load_data


# --------------------
# Core functions
# --------------------
def compute_flow(edges):
    """Compute total outflow (provider strength) and inflow (receiver strength)."""
    outflow = edges.groupby("provider")["score"].sum().sort_values(ascending=False)
    inflow = edges.groupby("receiver")["score"].sum().sort_values(ascending=False)
    return outflow, inflow


def select_top_nodes(outflow, inflow, top_n=10):
    """Select top providers and receivers."""
    return outflow.head(top_n).index, inflow.head(top_n).index


def filter_edges(edges, top_providers, top_receivers):
    """Filter edges to top providers and receivers."""
    return edges[
        edges["provider"].isin(top_providers) & edges["receiver"].isin(top_receivers)
    ]


def pivot_edges_to_matrix(edges):
    """Pivot edge list to provider-receiver matrix."""
    matrix = edges.pivot_table(
        index="provider", columns="receiver", values="score", fill_value=0
    )
    return matrix.astype(int)


# --------------------
# Heatmap
# --------------------
def plot_heatmap(matrix, title, output_prefix=None, figsize=(10, 6)):
    """Plot a heatmap from a pivoted matrix."""
    fig = plt.figure(figsize=figsize)
    sns.heatmap(matrix, cmap="viridis", annot=True, fmt="d")
    plt.title(title)
    plt.xlabel("Receiver")
    plt.ylabel("Provider")
    plt.tight_layout()

    # Save figure if output_prefix is provided
    if output_prefix is not None:
        plt.savefig(
            INTERACTION_ANALYSIS_PLOTS_DIR / f"{output_prefix}_heatmap.png", dpi=300
        )

    plt.close(fig)


# --------------------
# Process
# --------------------
def process_edge_file(edges_file, meta, output_prefix, heatmap_title):
    edges = load_data(edges_file, filetype="csv")

    # Compute flows
    outflow, inflow = compute_flow(edges)
    top_providers, top_receivers = select_top_nodes(outflow, inflow)

    # Filter edges & pivot
    sub_edges = filter_edges(edges, top_providers, top_receivers)
    matrix = pivot_edges_to_matrix(sub_edges)

    # Plot heatmap
    plot_heatmap(matrix, title=heatmap_title, output_prefix=output_prefix)

    # Metadata join & save
    top_providers_meta = pd.DataFrame({"patric_id": top_providers}).merge(
        meta, on="patric_id", how="left"
    )
    top_receivers_meta = pd.DataFrame({"patric_id": top_receivers}).merge(
        meta, on="patric_id", how="left"
    )

    top_providers_meta.to_csv(
        INTERACTION_ANALYSIS_PLOTS_DIR / f"{output_prefix}_top10_providers.tsv",
        sep="\t",
        index=False,
    )
    top_receivers_meta.to_csv(
        INTERACTION_ANALYSIS_PLOTS_DIR / f"{output_prefix}_top10_receivers.tsv",
        sep="\t",
        index=False,
    )

    print(f"Processed {edges_file}: heatmap + top providers/receivers saved")


# --------------------
# Main
# --------------------
def main():
    # Output directorty
    INTERACTION_ANALYSIS_PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    # Load metadata
    meta = load_data(METABOLIC_POTENTIAL_TSV, filetype="tsv")

    # Process each dataset
    process_edge_file(CYANO_HALO_EDGES_CSV, meta, "cyano_halo", "Heatmap: Cyano-Halo")
    process_edge_file(OUTLIERS_L2H_EDGES_CSV, meta, "low2high", "Heatmap: Low to High")
    process_edge_file(OUTLIERS_H2L_EDGES_CSV, meta, "high2low", "Heatmap: High to Low")


if __name__ == "__main__":
    main()
