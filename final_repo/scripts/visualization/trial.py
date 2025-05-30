import pandas as pd
import ast
import matplotlib.pyplot as plt

# Load data
data = pd.read_csv("/home/annavouk/master_thesis/final_repo/scripts/analysis/compounds_seed_nonseed_with_kegg.csv")

# Function to safely parse KEGG Module column (string of list)
def parse_kegg_modules(x):
    try:
        # If already a list, return as is
        if isinstance(x, list):
            return [m.strip() for m in x]
        # Parse string to list
        modules = ast.literal_eval(x)
        return [m.strip() for m in modules]
    except:
        return []

# Apply parsing
data['KEGG Module List'] = data['KEGG Module'].apply(parse_kegg_modules)

# Count KEGG modules per seed
modules_per_seed = data.groupby('ModelSEED ID')['KEGG Module List'].apply(lambda lists: set(sum(lists, []))).reset_index()
modules_per_seed['Num_KEGG_Modules'] = modules_per_seed['KEGG Module List'].apply(len)

print(modules_per_seed['Num_KEGG_Modules'])
print("Summary statistics of KEGG modules per seed:")
print(modules_per_seed['Num_KEGG_Modules'].describe())

plt.hist(modules_per_seed['Num_KEGG_Modules'], bins=range(1, modules_per_seed['Num_KEGG_Modules'].max()+2), edgecolor='black')
plt.xlabel('Number of KEGG Modules')
plt.ylabel('Number of Seeds')
plt.title('Distribution of KEGG Modules per Seed Compound')
plt.show()

# 2. To find which modules are most frequently complemented by seeds:
# Explode the KEGG Module list to one row per seed-module
exploded = data[['ModelSEED ID', 'KEGG Module List']].explode('KEGG Module List')

# Count how many unique seeds per module
module_seed_counts = exploded.groupby('KEGG Module List')['ModelSEED ID'].nunique().reset_index()
module_seed_counts.columns = ['KEGG Module', 'Num_Seed_Compounds']

print("Top 10 KEGG modules by number of seed compounds:")
print(module_seed_counts.sort_values('Num_Seed_Compounds', ascending=False).head(10))

# Plot top 10 KEGG modules
top10 = module_seed_counts.sort_values('Num_Seed_Compounds', ascending=False).head(10)
plt.bar(top10['KEGG Module'], top10['Num_Seed_Compounds'])
plt.xticks(rotation=45)
plt.ylabel('Number of Seed Compounds')
plt.title('Top 10 KEGG Modules by Seed Compound Count')
plt.show()
