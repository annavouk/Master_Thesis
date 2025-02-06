import pandas as pd

def get_seeds_subset(df, patric_ids):
    existing_ids = [pid for pid in patric_ids if pid in df.index]
    return df.loc[existing_ids]

#Example
if __name__ == "__main__":
    seeds_binary_df = pd.read_pickle('seeds_binary_per_patric.pckl')

    #Example list of PATRIC IDs to filter
    patric_ids_of_interest = ['100.11', '1000565.3', '1000000.0']  # Note: '1000000.0' does not exist

    subset_df = get_seeds_subset(seeds_binary_df, patric_ids_of_interest)

    print("\nSubset:")
    print(subset_df)
