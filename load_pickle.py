import pickle

# pickle_file = 'seeds_binary_per_patric.pckl'
pickle_file = 'seeds_per_genome_result_df.pkl'

with open(pickle_file, 'rb') as f:
	data = pickle.load(f)

print(data)
