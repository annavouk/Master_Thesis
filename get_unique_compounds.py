import pandas as pd

df1 = pd.read_pickle("seeds_binary_per_patric.pckl")
df2 = pd.read_pickle("non_seeds_binary_per_patric.pckl")

compound_columns_1 = [col for col in df1.columns if col.startswith('cpd')]
compound_columns_2 = [col for col in df2.columns if col.startswith('cpd')]

compounds = pd.concat([pd.Series(compound_columns_1), pd.Series(compound_columns_2)], ignore_index=True)
compounds = compounds.drop_duplicates()
compounds = compounds.sort_values()

compounds.to_csv("unique_compounds.csv", index=False, header=['Compound'])

print(compounds.head())
print(compounds.tail())
