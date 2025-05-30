import pandas as pd
import ast
from collections import defaultdict

# Load CSV
df = pd.read_csv("/home/annavouk/master_thesis/final_repo/scripts/analysis/compounds_seed_nonseed_with_kegg.csv")

# Parse string list of modules into actual lists
df['KEGG Module'] = df['KEGG Module'].apply(lambda x: [m.strip() for m in ast.literal_eval(x)])

# Prepare dictionary to accumulate seed genome counts per KEGG module
module_seed_genome_counts = defaultdict(int)

# Iterate through metabolites
for _, row in df.iterrows():
    seed_count = row['seed_genome_count']
    modules = row['KEGG Module']

    # Add seed_count to each module this metabolite participates in
    for module in modules:
        module_seed_genome_counts[module] += seed_count

# Convert to DataFrame and sort descending
module_distribution = pd.DataFrame(
    module_seed_genome_counts.items(),
    columns=['KEGG Module', 'Total Seed Genome Count']
).sort_values(by='Total Seed Genome Count', ascending=False)

print(module_distribution)
