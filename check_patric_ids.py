import csv

def read_single_column(filename):
    with open(filename, 'r') as file:
            return {line.strip() for line in file if line.strip()}  # Remove empty lines

def read_second_column(filename):
    with open(filename, 'r') as file:
            reader = csv.reader(file, delimiter='\t')
            # Create a dictionary with the ID as the key and both columns as values
            return {row[1].strip(): row[0].strip() for row in reader if len(row) > 1 and row[1].strip()}

def write_ids_to_file(filename, common_ids, data):
    with open(filename, 'w') as file:
            for id in sorted(common_ids):
                if id in data:
                    file.write(f"{id}\t{data[id]}\n")  # Write ID and corresponding second column

# File paths
file1 = 'patric_ids.csv'
file2 = 'gtdb2patricIds.tsv'
output_file = 'common_patric_ids.txt'

# Read data
ids1 = read_single_column(file1)  # First file (single column)
ids2 = read_second_column(file2)  # Second file (second column as key and first column as value)

# Find common IDs
common = ids1 & ids2.keys()  # Find common IDs between file1 and the second column of file2

# Save common IDs with the first and second columns from file2
write_ids_to_file(output_file, common, ids2)

print(f"Common IDs saved to {output_file}.")
