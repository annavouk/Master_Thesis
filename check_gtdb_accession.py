import pandas as pd

patric_file = "patric_ids_with_gtdb_taxonomy.csv"
gtdb_file = "gtdb_metadata_r207.tsv"

patric_df = pd.read_csv(patric_file)
gtdb_df = pd.read_csv(gtdb_file, sep="\t")

patric_assembly_accessions = set(patric_df["assembly_accession"])
gtdb_assembly_accessions = set(gtdb_df["accession"])

missing_accessions = patric_assembly_accessions - gtdb_assembly_accessions

if not missing_accessions:
    print("All assembly_accessions from patric_ids_with_gtdb_taxonomy.csv exist in gtdb_metadata_r207.tsv.")
else:
    print(f" {len(missing_accessions)} assembly_accessions are missing:")
    print("\n".join(missing_accessions))
