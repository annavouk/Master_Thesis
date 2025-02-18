import pandas as pd

df1 = pd.read_pickle("seeds_binary_per_patric.pckl")
df2 = pd.read_pickle("non_seeds_binary_per_patric.pckl")

compound_columns_df1 = [col for col in df1.columns if col.startswith('cpd')]
compound_columns_df2 = [col for col in df2.columns if col.startswith('cpd')]

compounds_df = pd.concat([pd.Series(compound_columns_df1), pd.Series(compound_columns_df2)], ignore_index=True)
compounds_df = compounds_df.drop_duplicates()

compounds_df.to_csv("unique_compounds.csv", index=False, header=False)

print(compounds_df.head())
print(compounds_df.tail())
