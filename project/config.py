from pathlib import Path

# Project Root
PROJECT_ROOT = Path(__file__).resolve().parent

# Directory Structure
INPUT_DIR = PROJECT_ROOT / "input"
OUTPUT_DIR = PROJECT_ROOT / "output"
PLOTS_DIR = OUTPUT_DIR / "plots"

METADATA_DIR = PROJECT_ROOT / "metadata"
METADATA_RAW_DIR = METADATA_DIR / "raw"

SCRIPTS_DIR = PROJECT_ROOT / "scripts"
UTILS_DIR = PROJECT_ROOT / "utils"

GTDB_TREES_DIR = PROJECT_ROOT / "gtdb_trees"

# Input Files
SEEDS_PICKLE = INPUT_DIR / "seeds_binary_per_patric.pckl"
NON_SEEDS_PICKLE = INPUT_DIR / "non_seeds_binary_per_patric.pckl"

# Metadata Files
GTDB_METADATA = METADATA_RAW_DIR / "gtdb_metadata_r207.tsv"
PATRIC_METADATA_JSON = METADATA_RAW_DIR / "patric_ids_metadata.json"
FAILED_PATRIC_IDS_TXT = METADATA_RAW_DIR / "failed_patric_ids.txt"
COMPACT_METADATA = METADATA_DIR / "compact_metadata.csv"


# Compound-level metadata
COMPOUNDS_TSV = METADATA_RAW_DIR / "compounds.tsv"
REACTIONS_TSV = METADATA_RAW_DIR / "reactions.tsv"
PATHWAYS_TSV = METADATA_RAW_DIR / "KEGG.pathways"
KEGG_MODULE_MAPPING = METADATA_RAW_DIR / "seedId_keggId_module.tsv"
COMPOUND_SUMMARY_TSV = METADATA_DIR / "compound_summary.tsv"
MODULE_MAP = METADATA_RAW_DIR/ "module_map_pairs.tsv"

# Output Files
METABOLIC_POTENTIAL_1 = OUTPUT_DIR / "metabolic_potential_summary.csv"
GENOME_PATHWAY_COVERAGE = OUTPUT_DIR / "genome_pathway_coverage.tsv"
AMINOACIDS = OUTPUT_DIR / "amino_acid_compounds.tsv"
