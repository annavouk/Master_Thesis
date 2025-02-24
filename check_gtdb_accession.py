import pandas as pd

patric_file = "patric_ids_with_gtdb_taxonomy.csv" # data from GTDB2 file (every PATRIC ID of interest with accession number) and GTDB_accession_taxonomy.csv
gtdb_file = "gtdb_accession_taxonomy.csv" # metadata from GTDB (gtdb_metadata_r207.tsv)
output_file = "missing_accessions.txt"

patric_df = pd.read_csv(patric_file, low_memory=False)
gtdb_df = pd.read_csv(gtdb_file, low_memory=False)

patric_assembly_accessions = set(patric_df["assembly_accession"])
gtdb_assembly_accessions = set(gtdb_df["assembly_accession"])

missing_accessions = patric_assembly_accessions - gtdb_assembly_accessions

if not missing_accessions:
    print("All assembly_accessions from patric_ids_with_gtdb_taxonomy.csv exist in gtdb_metadata_r207.tsv.")
else:
    print(f"{len(missing_accessions)} assembly_accessions are missing. Saved to {output_file}")

    with open(output_file, "w") as f:
        f.write("\n".join(missing_accessions))
