import pandas as pd

def seeds_per_genome(df):
	df['Total_Seeds'] = df.apply(lambda row: (row == 1).sum(), axis=1)
	return df

if __name__ == "__main__":

	df = pd.read_pickle('seeds_binary_per_patric.pckl')

	result_df = seeds_per_genome(df)

	print(result_df)
