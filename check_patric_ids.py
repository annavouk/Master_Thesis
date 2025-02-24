import csv

def read_single_column(filename):
    with open(filename, 'r') as file:
        return {line.strip() for line in file if line.strip()}  # Remove empty lines

def read_second_column(filename):
    with open(filename, 'r') as file:
        reader = csv.reader(file, delimiter='\t')  # Tab-separated file
        return {row[1].strip() for row in reader if len(row) > 1 and row[1].strip()}  # Extract second column

def write_ids_to_file(filename, ids):
    with open(filename, 'w') as file:
        for id in sorted(ids):
            file.write(f"{id}\n")

# File paths
file1 = 'patric_ids.csv'
file2 = 'gtdb2patricIds.tsv'
output_file = 'differences.txt'

# Read data
ids1 = read_single_column(file1)  # First file (single column)
ids2 = read_second_column(file2)  # Second file (second column)

# Find differences
differences = ids1 - ids2  # Keep only IDs in patric_ids.csv but not in gtdb2patricIds.tsv

# Save differences
write_ids_to_file(output_file, differences)

print(f"Differences saved to {output_file}.")
