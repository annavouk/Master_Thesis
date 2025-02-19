import pandas as pd

cpds = pd.read_csv("unique_compounds.csv")
cpds_list = cpds["Compound"].astype(str).tolist()
reactions = pd.read_csv("reactions_compact.tsv", sep='\t', low_memory=False).astype(str)
compounds_reactions = {}

for cpd in cpds_list:
	for row in reactions.itertuples(index=False, name=None):
		if any(cpd in str(x) for x in row):
			compounds_reactions.setdefault(cpd, []).append(list(row))


df = pd.DataFrame([(k, v) for k, v in compounds_reactions.items()], columns=["Compound", "Reactions"])

df.to_csv("compounds_reactions.tsv", sep='\t', index=False)

print(df.head())
print(df.tail())
