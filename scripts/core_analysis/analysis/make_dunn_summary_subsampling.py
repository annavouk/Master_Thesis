import pandas as pd
from io import StringIO

pairs = []
with open("output/dunn_posthoc_subsampling.tsv") as f:
    rep = None
    metric = None
    block = []
    for line in f:
        if line.startswith("# Rep"):
            # πχ "# Rep 0 Metric Total_Seeds"
            parts = line.strip().split()
            rep = int(parts[2])
            metric = parts[-1]
            block = []
        elif line.strip():
            block.append(line)
        elif block:
            df = pd.read_csv(StringIO("".join(block)), sep="\t", index_col=0)
            melted = df.stack().reset_index()
            melted.columns = ["Group1", "Group2", "p_value"]
            melted = melted[melted["Group1"] != melted["Group2"]]
            melted["Pair"] = melted.apply(
                lambda r: " vs ".join(sorted([r["Group1"], r["Group2"]])), axis=1
            )
            melted["Metric"] = metric
            melted["Rep"] = rep
            pairs.append(melted[["Rep", "Metric", "Pair", "p_value"]])
            block = []

# Combine reps
dunn_subs_long = pd.concat(pairs)

# Summary
summary = (
    dunn_subs_long
    .groupby(["Metric", "Pair"])
    .agg(
        Significant_reps=("p_value", lambda x: f"{(x<0.05).sum()} / {len(x)}"),
        Percent_significant=("p_value", lambda x: 100*(x<0.05).mean()),
        Median_p_sig=("p_value", lambda x: x[x<0.05].median() if any(x<0.05) else None),
        Min_p_sig=("p_value", lambda x: x[x<0.05].min() if any(x<0.05) else None),
        Max_p_sig=("p_value", lambda x: x[x<0.05].max() if any(x<0.05) else None),
    )
    .reset_index()
)

summary.to_csv("output/dunn_posthoc_summary.tsv", sep="\t", index=False)
print("Saved Dunn’s subsampling summary")
