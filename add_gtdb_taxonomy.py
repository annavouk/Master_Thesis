import csv

def read_common_patric_ids(filename):
	with open(filename, 'r') as file:
		return [tuple(line.strip().split("\t")) for line in file if line.strip()]

def read_gtdb_metadata(filename):
	gtdb_dict = {}
	with open(filename, 'r') as file:
		reader = csv.reader(file)
		next(reader)

		for row in reader:
			if len(row) >= 3:
				accession1 = row[0].strip()
				accession2 = row[2].strip()
				taxonomy = row[1].strip()

				gtdb_dict[accession1] = taxonomy
				gtdb_dict[accession2] = taxonomy
	return gtdb_dict

def write_output(filename, patric_data, gtdb_dict):
	with open(filename, 'w') as file:
		for patric_id, accession in patric_data:
			taxonomy = gtdb_dict.get(accession, "Not Found")
			file.write(f"{patric_id}\t{accession}\t{taxonomy}\n")

common_patric_file = 'common_patric_ids.txt'
gtdb_metadata_file = 'gtdb_metadata_compact.csv'
output_file = 'common_patric_ids_with_taxonomy.txt'

patric_data = read_common_patric_ids(common_patric_file)
gtdb_data = read_gtdb_metadata(gtdb_metadata_file)

missing = [accession for _, accession in patric_data if accession not in gtdb_data]
print(f"Missing {len(missing)} accessions:", missing[:10])

write_output(output_file, patric_data, gtdb_data)

print(f"{output_file} has {len(patric_data)} rows.")

