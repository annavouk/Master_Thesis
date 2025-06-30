import pandas as pd

df = pd.read_csv('metadata/compound_summary.tsv', sep='\t', low_memory=False)

amino_acids = [
    'alanine', 'arginine', 'asparagine', 'aspartate', 'aspartic acid', 'cysteine',
    'glutamine', 'glutamate', 'glutamic acid', 'glycine', 'histidine', 'isoleucine',
    'leucine', 'lysine', 'methionine', 'phenylalanine', 'proline', 'serine',
    'threonine', 'tryptophan', 'tyrosine', 'valine'
]

regex_pattern = '|'.join(amino_acids)

aa_df = df[df['compound_name'].str.contains(regex_pattern, case=False, na=False)]

print(aa_df.head())

aa_df.to_csv('output/amino_acid_compounds.tsv', sep='\t', index=False)



