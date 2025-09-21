"""
Metabolic Interaction Analysis

This script builds feeding matrices (provider-receiver) for specific genome subsets,
based on metadata filters (biome, phylum, ratio thresholds, custom IDs).

Modes:
- In-memory (default): computes full provider-receiver matrices in RAM, exports edge lists, summary.
- Blockwise (--blockwise): processes large datasets in blocks, saves matrices + edge lists directly to CSVs.

Usage:
    python3 scripts/core_analysis/analysis/provider_receiver.py --subset 200 
    python3 scripts/core_analysis/analysis/provider_receiver.py --blockwise --block-size 1000 --out-prefix all
    python3 provider_receiver.py --by-biome Soil
    python3 provider_receiver.py --taxon-rank phylum --taxa Halobacteriota Cyanobacteria --cross-taxa-only
    python3 provider_receiver.py --taxon-rank family --taxa Planctomycetaceae Verrucomicrobiaceae

Note:
    Blockwise mode processes the entire dataset in chunks (ignores downstream filtering by taxonomy/biome).
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

import argparse
import pandas as pd
import numpy as np

from config import(
    SEEDS_PICKLE,
    NON_SEEDS_PICKLE,
    COMPACT_METADATA_WITH_BIOME_TSV,
    METABOLIC_POTENTIAL_TSV,
    OUTPUT_DIR,
    )
from utils import load_data, parse_taxonomy


# ------------------------
# Helpers for output paths
# ------------------------
def make_output_paths(prefix="feeding", blockwise=False, subset=None):
    """Return dict of output paths with suffixes depending on mode."""
    suffix = []
    if blockwise:
        suffix.append("blockwise")
    if subset:
        suffix.append(f"subset{subset}")
    suf = "_" + "_".join(suffix) if suffix else ""

    return {
        "matrix": OUTPUT_DIR / f"{prefix}_matrix{suf}.csv",
        "edges": OUTPUT_DIR / f"{prefix}_edgelist{suf}.csv",
    }


# ------------------------
# Core computations
# ------------------------
def compute_feeding_matrix(seeds_df, non_seeds_df):
    """Return feeding matrices (in-memory)."""
    common_cpds = seeds_df.columns.intersection(non_seeds_df.columns)     # compounds that are both seeds and non-seeds
    seeds_bin = seeds_df[common_cpds].values
    nonseeds_bin = non_seeds_df[common_cpds].values

    abs_mat = np.dot(nonseeds_bin, seeds_bin.T)    # raw counts of shared compounds

    return abs_mat


def feeding_matrix_to_edgelist(feeding_df, threshold):
    """Convert feeding matrix (DataFrame) to edge list with scores >= threshold."""
    feeding_df = feeding_df.rename_axis(index="provider", columns="receiver").copy()

    edges = feeding_df.stack().reset_index(name="score")
    edges = edges[(edges["score"] >= threshold) & (edges["provider"] != edges["receiver"])]

    return edges


def compute_feeding_matrix_blockwise(seeds_df, non_seeds_df, paths, raw_thr=1, block_size=1000):
    """Compute feeding matrix block-by-block (for interactions across the full dataset)."""
    common_cpds = seeds_df.columns.intersection(non_seeds_df.columns)
    seeds_bin = seeds_df[common_cpds].values.astype(np.uint8)
    nonseeds_bin = non_seeds_df[common_cpds].values.astype(np.uint8)

    n_providers = nonseeds_bin.shape[0]

    # Init CSVs
    pd.DataFrame(columns=seeds_df.index).to_csv(paths["matrix"], index=False)
    pd.DataFrame(columns=["provider", "receiver", "score"]).to_csv(paths["edges"], index=False)

    # Process in chunks
    for start in range(0, n_providers, block_size):
        end = min(start + block_size, n_providers)
        block = nonseeds_bin[start:end, :]

        abs_block = np.dot(block, seeds_bin.T)

        providers = non_seeds_df.index[start:end]
        receivers = seeds_df.index

        # Save matrices
        feeding_matrix_df = pd.DataFrame(abs_block, index=providers, columns=receivers)
        feeding_matrix_df.index.name = None
        feeding_matrix_df.columns.name = None
        feeding_matrix_df.to_csv(paths["matrix"], mode="a", header=False)

        # Save edges
        edges = (
            feeding_matrix_df
            .rename_axis(index="provider", columns="receiver")
            .stack()
            .reset_index(name="score")
            .query("score >= @raw_thr and provider != receiver")
            )

        edges.to_csv(paths["edges"], mode="a", header=False, index=False)
        
        print(f"Processed block {start}:{end} / {n_providers}")

    print(f"\nBlockwise outputs saved:")
    for k, v in paths.items():
        print(f" - {v}")


# ------------------------
# Filtering helpers
# ------------------------
def filter_cross_taxa_edges(edges, df, rank):
    """Keep only edges where provider and receiver belong to different taxa at given rank."""
    tax_map = df.set_index("patric_id")[rank].to_dict()
    edges[f"{rank}_provider"] = edges["provider"].map(tax_map)
    edges[f"{rank}_receiver"] = edges["receiver"].map(tax_map)
    return edges[edges[f"{rank}_provider"] != edges[f"{rank}_receiver"]]


# ------------------------
# Run modes
# ------------------------
def run_in_memory(seeds_df, non_seeds_df, df, args, paths):
    abs_mat = compute_feeding_matrix(seeds_df, non_seeds_df)
    feeding_raw = pd.DataFrame(abs_mat, index=non_seeds_df.index, columns=seeds_df.index)

    feeding_raw.to_csv(paths["matrix"])
    print("Feeding matrices saved (in-memory mode).")

    edges = feeding_matrix_to_edgelist(feeding_raw, threshold=args.raw_threshold)

    if args.cross_taxa_only:
        edges = filter_cross_taxa_edges(edges, df, args.taxon_rank)

    edges.to_csv(paths["edges"], index=False)
    print("Edge lists exported.")


def run_blockwise(seeds_df, non_seeds_df, args, paths):
    compute_feeding_matrix_blockwise(
        seeds_df, non_seeds_df,
        paths=paths,
        raw_thr=args.raw_threshold,
        block_size=args.block_size
    )


# ------------------------
# Debug
# ------------------------
def debug_report(df, seeds_df, non_seeds_df, label=""):
    print(f"\n[DEBUG {label}]")
    print(" - genomes left in metadata:", df.shape[0])
    print(" - seeds_df shape:", seeds_df.shape)
    print(" - non_seeds_df shape:", non_seeds_df.shape)
    print(" - common compounds:", len(seeds_df.columns.intersection(non_seeds_df.columns)))
    if "main_biome" in df.columns:
        print(" - biomes:", df["main_biome"].unique())
    if "phylum" in df.columns:
        print(" - phyla:", df["phylum"].unique()[:10])


# ------------------------
# Main
# ------------------------
def main():
    parser = argparse.ArgumentParser(description="Metabolic Interaction Analysis (Feeding Matrix, raw only)")
    parser.add_argument("--low-high-ratio", nargs=2, type=float, default=None,
                        help="Compare genomes from low vs high ratio groups (e.g. --low-high-ratio 0.2 0.4)")
    parser.add_argument("--similar-ratio", nargs=2, type=float, default=None,
                        help="Restrict to genomes with ratio between two thresholds (e.g. --similar-ratio 0.30 0.32)")
    parser.add_argument("--taxon-rank", type=str, default="phylum",
                        help="Taxonomic rank to filter (e.g. phylum, class, order, family, genus)")
    parser.add_argument("--taxa", nargs="+", type=str, default=None,
                        help="List of taxa to include (e.g. --taxa Halobacteriota Cyanobacteria)")
    parser.add_argument("--cross-taxa-only", action="store_true",
                        help="Keep only edges between different taxa at the chosen rank")
    parser.add_argument("--by-biome", type=str, default=None,
                        help="Restrict analysis to genomes from this biome (e.g. Soil, Marine, Freshwater)")
    parser.add_argument("--debug", action="store_true",
                        help="Print debug info after filtering")
    parser.add_argument("--raw-threshold", type=int, default=1,
                        help="Threshold for feeding edge list (default=1)")
    parser.add_argument("--subset", type=int, default=None, help="Subset of genomes for testing")
    parser.add_argument("--blockwise", action="store_true", help="Use blockwise computation")
    parser.add_argument("--block-size", type=int, default=1000, help="Block size for blockwise mode")
    parser.add_argument("--out-prefix", type=str, default=None,
                        help="Prefix for output files (if not given, auto-generated)")

    args = parser.parse_args()

    # Load Pickle dataframes
    seeds_df = load_data(SEEDS_PICKLE, filetype="pickle")
    non_seeds_df = load_data(NON_SEEDS_PICKLE, filetype="pickle")

    # Load metadata
    biomes_df = load_data(COMPACT_METADATA_WITH_BIOME_TSV, filetype="tsv")
    ratio_df = load_data(METABOLIC_POTENTIAL_TSV, filetype="tsv")

    # Keep only "main_biome" from metadata
    biomes_df = biomes_df[["patric_id", "main_biome"]]

    # Merge + parse taxonomy
    df = ratio_df.merge(biomes_df, on="patric_id", how="left")
    df = parse_taxonomy(df)

    # Apply filters on metadata
    if args.by_biome:
        df = df[df["main_biome"] == args.by_biome]

    if args.taxa:
        df = df[df[args.taxon_rank].isin(args.taxa)]

    if args.low_high_ratio:
        low_thr, high_thr = args.low_high_ratio
        df = df[(df["Ratio"] <= low_thr) | (df["Ratio"] >= high_thr)]

    if args.similar_ratio:
        low_thr, high_thr = args.similar_ratio
        df = df[(df["Ratio"] >= low_thr) & (df["Ratio"] <= high_thr)]

    if args.subset:
        df = df.head(args.subset)

    # Subset matrices
    selected_ids = df["patric_id"].astype(str).tolist()
    seeds_df = seeds_df.loc[seeds_df.index.isin(selected_ids)]
    non_seeds_df = non_seeds_df.loc[non_seeds_df.index.isin(selected_ids)]

    if args.debug:
        debug_report(df, seeds_df, non_seeds_df, "after filters")

    # Build output prefix
    if args.out_prefix:
        prefix = args.out_prefix
    else:
        suffix_parts = []
        if args.low_high_ratio:
            low_thr, high_thr = args.low_high_ratio
            suffix_parts.append(f"ratio{low_thr}-{high_thr}")
        if args.similar_ratio:
            low_thr, high_thr = args.similar_ratio
            suffix_parts.append(f"similar{low_thr}-{high_thr}")
        if args.taxa:
            taxa_str = "-".join(args.taxa)
            suffix_parts.append(f"{args.taxon_rank}-{taxa_str}")
        if args.cross_taxa_only:
            suffix_parts.append("cross")
        if args.by_biome:
            suffix_parts.append(f"biome-{args.by_biome}")

        prefix = "feeding" + ("_" + "_".join(suffix_parts) if suffix_parts else "")

    paths = make_output_paths(prefix=prefix, blockwise=args.blockwise, subset=args.subset)

    # Run computation
    if args.blockwise:
        run_blockwise(seeds_df, non_seeds_df, args, paths)
    else:
        run_in_memory(seeds_df, non_seeds_df, df, args, paths)


if __name__ == "__main__":
    main()
