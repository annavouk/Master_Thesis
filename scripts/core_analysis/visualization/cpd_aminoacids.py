"""
Analyze amino acid - Auxotrophy & Prototrophy profiles across genomes

Generate plots:
    - Barplot of auxotrophy/prototrophy frequency per amino acid
    - Histogram of auxotrophies/prototrophies per genome
    - Stacked barplots of auxotrophy patterns (%) per taxonomic rank
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[3]))

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm

from config import (
    SEEDS_PICKLE,  # input
    NON_SEEDS_PICKLE,  # input
    COMPACT_METADATA_TSV,  # metadata
    AMINOACIDS_TSV,  # aa metadata
    OUTPUT_DIR,  # output dir
    AMINOACIDS_PLOTS_DIR,  # plots dir
)
from utils import load_data, split_and_clean_taxonomy


# ------------------------
# Utilities
# ------------------------
def ensure_str_index_and_columns(df):
    """Make index and columns string-typed for safe alignment/merges."""
    df = df.copy()
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def ensure_binary(df):
    """Clip values to {0,1} and use compact dtype."""
    return df.fillna(0).clip(0, 1).astype("int8")


# ------------------------
# Auxotrophy analysis
# ------------------------
def analyze_auxotrophy(seed_matrix, nonseed_matrix, aa_df, metadata_df):
    """Compute auxotrophy/prototrophy distribution per genome with taxonomic aggregation."""
    # Normalize typing
    seed_matrix = ensure_str_index_and_columns(seed_matrix)
    nonseed_matrix = ensure_str_index_and_columns(nonseed_matrix)

    # Reference list of amino acids (ModelSEED IDs) from AMINOACIDS_TSV (extracted from KEGG)
    aa_cpd_ids = aa_df["ModelSEED_ID"].astype(str).tolist()

    # Union (present in both seed and non-seed matrices)
    aa_in_seed = [c for c in aa_cpd_ids if c in seed_matrix.columns]
    aa_in_non = [c for c in aa_cpd_ids if c in nonseed_matrix.columns]
    common_aa = sorted(set(aa_in_seed).union(aa_in_non))

    # Subset and enforce binary
    seed_aa_matrix = ensure_binary(seed_matrix.reindex(columns=common_aa, fill_value=0))
    nonseed_aa_matrix = ensure_binary(
        nonseed_matrix.reindex(columns=common_aa, fill_value=0)
    )

    # Genome-level auxotrophy/prototrophy counts
    auxotrophy_per_genome = seed_aa_matrix.sum(axis=1)
    prototrophy_per_genome = nonseed_aa_matrix.sum(axis=1)

    summary_df = pd.DataFrame(
        {"Auxotrophy": auxotrophy_per_genome, "Prototrophy": prototrophy_per_genome}
    )
    summary_df.index.name = "PATRIC"

    # Amino acid-level auxotrophy/prototrophy counts
    seed_counts = seed_aa_matrix.sum(axis=0)
    nonseed_counts = nonseed_aa_matrix.sum(axis=0)

    counts_df = pd.DataFrame(
        {
            "Auxotrophy": seed_counts,
            "Prototrophy": nonseed_counts,
        }
    )

    # Replace ModelSEED IDs with compound names
    labels = aa_df.set_index("ModelSEED_ID")["compound_name"].astype(str)
    names = labels.reindex(counts_df.index)
    # If compound name is missing, use ModelSEED ID as fallback
    fallback = pd.Series(counts_df.index.astype(str), index=counts_df.index)
    names = names.fillna(fallback)

    counts_df.index = names
    counts_df.index.name = "Amino Acid"

    # Add prefix to amino acid names
    seed_pref = seed_aa_matrix.add_prefix("seed__")
    nonseed_pref = nonseed_aa_matrix.add_prefix("nonseed__")
    both = pd.concat([seed_pref, nonseed_pref], axis=1)
    both.index = both.index.str.strip()
    metadata_df["patric_id"] = metadata_df["patric_id"]

    # Merge with taxonomy
    merged = both.reset_index(names="PATRIC").merge(
        metadata_df[["patric_id", "gtdb_taxonomy"]],
        left_on="PATRIC",
        right_on="patric_id",
        how="left",
    )
    taxonomy_df = split_and_clean_taxonomy(merged, "gtdb_taxonomy")
    merged = pd.concat([merged.drop(columns=["patric_id"]), taxonomy_df], axis=1)
    merged.set_index("PATRIC", inplace=True)

    return (
        summary_df,
        counts_df,
        merged,
        list(seed_pref.columns),
        list(nonseed_pref.columns),
    )


# ------------------------
# Taxonomic aggregation
# ------------------------
def aggregate_by_taxon(
    merged, seed_cols_pref, nonseed_cols_pref, ranks=("phylum", "family", "genus")
):
    """Compute auxotrophy/prototrophy counts and percentages per amino acid per taxonomic rank."""
    results = {}  # empty dictionary to store results

    for rank in ranks:
        # Group by taxonomic rank
        rank_group = merged.groupby(rank, dropna=False)

        # Count genomes per taxon
        aux_counts = rank_group[seed_cols_pref].sum(
            min_count=1
        )  # number of seed AAs per rank
        proto_counts = rank_group[nonseed_cols_pref].sum(min_count=1)

        rank_counts = merged.groupby(
            rank, dropna=False
        ).size()  # Number of total genomes per rank
        rank_counts = rank_counts.reindex(aux_counts.index).astype(float)

        # Calculate percentages per taxon to avoid bias from uneven sampling
        pct_aux = aux_counts.div(rank_counts, axis=0) * 100
        pct_proto = proto_counts.div(rank_counts, axis=0) * 100

        # Store results
        results[rank] = {
            "aux": {"counts": aux_counts, "percentages": pct_aux},
            "proto": {"counts": proto_counts, "percentages": pct_proto},
            "rank_counts": rank_counts,
        }

    return results


# ------------------------
# Visualization
# ------------------------
def plot_amino_acid_counts(counts_df, log_scale=False, save_path=None):
    """Barplot of auxotrophy vs prototrophy frequency per amino acid."""
    order = counts_df.sort_values("Auxotrophy", ascending=False).index
    fig, ax = plt.subplots(figsize=(10, 6))
    counts_df.loc[order].plot(kind="bar", ax=ax, edgecolor="black", alpha=0.9)
    plt.title("Auxotrophy/Prototrophy Frequency per Amino Acid")
    plt.ylabel("Number of Genomes")
    plt.xlabel("Amino Acid")

    if log_scale:
        ax.set_yscale("log")
        ax.set_ylabel("Number of Genomes (log scale)")

    plt.xticks(rotation=45, ha="right")

    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.2, -0.5),
        ncol=2,
        fontsize=9,
        frameon=False,
    )

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300)
        print(f"Saved plot: {save_path}")
    plt.close(fig)
    return None


def plot_summary_hist(summary_df, n_amino_acids, log_scale=False, save_path=None):
    """Histogram of auxotrophies and prototrophies per genome."""
    auxo = summary_df["Auxotrophy"]
    proto = summary_df["Prototrophy"]
    max_val = int(max(auxo.max(), proto.max(), n_amino_acids))
    bins = np.arange(0, max_val + 1, 1)  # 0–19

    fig = plt.figure(figsize=(10, 6))
    plt.hist(auxo, bins=bins, alpha=0.55, label="Auxotrophy", edgecolor="black")
    plt.hist(proto, bins=bins, alpha=0.55, label="Prototrophy", edgecolor="black")
    plt.xlabel("Number of Amino Acids")
    plt.ylabel("Number of Genomes")
    plt.title("Distribution of Auxotrophy and Prototrophy per Genome")

    if log_scale:
        plt.yscale("log")
        plt.ylabel("Number of Genomes (log scale)")
    plt.legend()
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300)
        print(f"Saved plot: {save_path}")
    plt.close(fig)
    return None


def plot_stacked(
    rank_summary,
    rank,
    aa_df,
    mode="percent",  # "percent" (per-taxon %) or "composition" (per-AA 100%)
    top_n=20,
    rank_counts=None,
    y_label=None,
    save_path=None,
):
    """Stacked barplot of auxotrophy per amino acid for the top-N taxa."""
    data = rank_summary.copy()

    # Select top-N taxa
    if rank_counts is not None:
        ordered_taxa = rank_counts.sort_values(ascending=False).index
        ordered_taxa = [t for t in ordered_taxa if t in data.index]  # Alignment
    else:
        ordered_taxa = data.sum(axis=1).sort_values(ascending=False).index
    top_taxa = list(ordered_taxa[:top_n])
    data = data.loc[top_taxa]

    # Replace prefixed column names with compound names
    cols = data.columns.astype(str)
    core_ids = [c.split("__", 1)[-1] if "__" in c else c for c in cols]

    labels = aa_df.assign(ModelSEED_ID=aa_df["ModelSEED_ID"].astype(str)).set_index(
        "ModelSEED_ID"
    )["compound_name"]
    new_cols = [labels.get(cid, cid) for cid in core_ids]
    data.columns = new_cols

    # Handle potential duplicate names after replacement
    if len(set(new_cols)) != len(new_cols):
        data = data.T.groupby(level=0).sum().T

    # Composition mode to ensure 100% per aminoacid
    if mode == "composition":
        col_sums = data.sum(axis=0).replace(0, np.nan)
        data = (data / col_sums) * 100
        data = data.fillna(0)

    # Colour for each taxon
    n = len(data.index)
    cmap = plt.get_cmap("tab20", n)
    colors = [cmap(i) for i in range(n)]

    # Plot (columns=AA, rows=taxa)
    ax = data.T.plot(
        kind="bar",
        stacked=True,
        figsize=(14, 9),
        color=colors,
    )

    fig = ax.figure

    if y_label is None:
        if mode == "composition":
            y_label = "Relative contribution per AA (%)"
        else:
            y_label = "Percentage of Genomes"

    ax.set_ylabel(y_label)
    title_mode = (
        "Composition (100% per AA)" if mode == "composition" else "Percent within taxon"
    )
    ax.set_title(f"Auxotrophy — {title_mode}\nTop {len(top_taxa)} {rank.capitalize()}")
    plt.xticks(rotation=45, ha="right")

    # Legend below
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.25),
        ncol=min(6, len(top_taxa)),
        fontsize=8,
        title=rank.capitalize(),
        frameon=False,
    )

    plt.tight_layout(rect=[0, 0.07, 1, 1])
    if save_path:
        fig.savefig(save_path, dpi=300)
        print(f"Saved plot: {save_path}")
    plt.close(fig)


# ------------------------
# Main
# ------------------------
def main():
    # Load data
    seeds = load_data(SEEDS_PICKLE, filetype="pickle")
    non_seeds = load_data(NON_SEEDS_PICKLE, filetype="pickle")
    aa_df = load_data(AMINOACIDS_TSV, filetype="tsv")
    metadata_df = load_data(COMPACT_METADATA_TSV, filetype="tsv", dtype=str)

    # Directory to save the plots
    plot_dir = Path(AMINOACIDS_PLOTS_DIR)
    plot_dir.mkdir(parents=True, exist_ok=True)

    # Analyze auxotrophy
    summary_df, counts_df, merged, seed_cols_pref, nonseed_cols_pref = (
        analyze_auxotrophy(seeds, non_seeds, aa_df, metadata_df)
    )

    # Taxonomic aggregation
    results = aggregate_by_taxon(merged, seed_cols_pref, nonseed_cols_pref)

    # Column for Valine prototrophy
    col_valine_proto = "nonseed__cpd00156"
    # Column for Isoleucine prototrophy
    col_isoleucine_proto = "nonseed__cpd00322"

    # Genomes prototrophic for Valine
    valine_producers = merged[merged[col_valine_proto] == 1]
    print("Genomes that produce Valine:", len(valine_producers))

    # Genomes prototrophic for Isoleucine
    isoleucine_producers = merged[merged[col_isoleucine_proto] == 1]
    print("Genomes that produce Isoleucine:", len(isoleucine_producers))

    # Per phylum count
    print(valine_producers["phylum"].value_counts().head(15))
    print(isoleucine_producers["phylum"].value_counts().head(15))

    # Save outputs
    counts_df.to_csv(
        OUTPUT_DIR / "amino_acid_auxotrophy_prototrophy_counts.tsv", sep="\t"
    )
    summary_df.to_csv(OUTPUT_DIR / "amino_acid_auxo_proto_per_genome.tsv", sep="\t")

    for rank, res in results.items():
        # counts + percentages
        aux_df = pd.concat(
            [res["aux"]["counts"], res["aux"]["percentages"]],
            axis=1,
            keys=["counts", "percentages"],
        )
        proto_df = pd.concat(
            [res["proto"]["counts"], res["proto"]["percentages"]],
            axis=1,
            keys=["counts", "percentages"],
        )

        combined = pd.concat(
            [aux_df, proto_df], axis=1, keys=["Auxotrophy", "Prototrophy"]
        )

        combined.columns = [
            f"{col.split('__')[-1]}_{lvl0}_{lvl1}"  # e.g. cpd00023_Auxotrophy_counts
            for lvl0, lvl1, col in combined.columns.to_flat_index()
        ]
        out_path = OUTPUT_DIR / f"auxo_proto_counts_pct_per_{rank}.tsv"
        combined.to_csv(out_path, sep="\t")
        print(f"Saved: {out_path}")

    # Plots
    plot_amino_acid_counts(
        counts_df, log_scale=True, save_path=plot_dir / "aa_counts.png"
    )
    n_aa = counts_df.shape[0]
    plot_summary_hist(
        summary_df,
        n_amino_acids=n_aa,
        log_scale=True,
        save_path=plot_dir / "aa_summary_hist.png",
    )

    for rank, res in results.items():
        order = res.get("rank_counts")
        plot_stacked(
            res["aux"]["counts"],
            rank,
            aa_df,
            mode="composition",
            top_n=20,
            rank_counts=order,
            save_path=plot_dir / f"aa_auxotrophy_COMPOSITION_{rank}.png",
        )


if __name__ == "__main__":
    main()
