import json
import csv

with open('patric_ids_metadata.json', 'r') as f:
	data = json.load(f)

patric_to_assembly = []

for patric_id, metadata in data.items():
	assembly_accession = metadata.get('assembly_accession', 'N/A')
	genome_length = metadata.get('genome_length', 'N/A')
	patric_to_assembly.append([patric_id, assembly_accession, genome_length])

output_file = 'patric_to_assembly.csv'

with open(output_file, 'w', newline='') as f:
	writer = csv.writer(f)
	writer.writerow(['patric_id', 'assembly_accession', 'genome_size'])
	writer.writerows(patric_to_assembly)

print(f"Data has been saved to {output_file}")
