import pandas as pd
import os

# Config
INPUT = "output/kruskal_subsampling.tsv"     # subsamples
FULL = "output/kruskal_biomes.tsv"   # full dataset
OUTPUT = "output/kruskal_summary.tsv"

# Load subsampling results
subs = pd.read_csv(INPUT, sep="\t")

# Load full if exists
if os.path.exists(FULL):
    full = pd.read_csv(FULL, sep="\t").reset_index().rename(columns={"index": "Metric"})
else:
    print(f"[WARN] Full results file not found: {FULL}")
    full = pd.DataFrame(columns=["Metric", "H_stat", "p_value"])

rows = []

for metric, g in subs.groupby("Metric"):
    sig = g[g["p"] < 0.05]
    n_sig = len(sig)
    perc_sig = (n_sig / len(g)) * 100 if len(g) > 0 else 0

    # Try to fetch from full dataset
    if metric in full["Metric"].values:
        H_full = float(full.loc[full["Metric"] == metric, "H_stat"].values[0])
        p_full = float(full.loc[full["Metric"] == metric, "p_value"].values[0])
    else:
        H_full, p_full = None, None

    rows.append({
        "Metric": metric,
        "H_stat_full": H_full,
        "p_value_full": p_full,
        "Significant_reps": f"{n_sig} / {len(g)}",
        "Percent_significant": perc_sig,
        "Median_H_sig": sig["H"].median() if n_sig > 0 else None,
        "Min_H_sig": sig["H"].min() if n_sig > 0 else None,
        "Max_H_sig": sig["H"].max() if n_sig > 0 else None,
        "Median_p_sig": sig["p"].median() if n_sig > 0 else None,
        "Min_p_sig": sig["p"].min() if n_sig > 0 else None,
        "Max_p_sig": sig["p"].max() if n_sig > 0 else None
    })

summary = pd.DataFrame(rows)

# Save
summary.to_csv(OUTPUT, sep="\t", index=False)
print(f"Saved Kruskal summary to {OUTPUT}")
print(summary)
