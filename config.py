from pathlib import Path


# Project Root
PROJECT_ROOT = Path(__file__).resolve().parent


# Directory Structure
INPUT_DIR = PROJECT_ROOT / "input"

OUTPUT_DIR = PROJECT_ROOT / "output"
PLOTS_DIR = OUTPUT_DIR / "plots"
EXPLORATORY_PLOTS_DIR = PLOTS_DIR / "exploratory_data_overview"
METABOLIC_POTENTIAL_PLOTS_DIR = PLOTS_DIR / "metabolic_potential"
BIOMES_PLOTS_DIR = PLOTS_DIR / "biomes_analysis"
SEED_NODE_FREQUENCY_PLOTS_DIR = PLOTS_DIR / "seed_node_frequency"

MAPS_DIR = PROJECT_ROOT / "maps"

METADATA_DIR = PROJECT_ROOT / "metadata"
METADATA_RAW_DIR = METADATA_DIR / "raw"

SCRIPTS_DIR = PROJECT_ROOT / "scripts"

UTILS_DIR = PROJECT_ROOT / "utils"

GTDB_TREES_DIR = PROJECT_ROOT / "gtdb_trees"

TEST_INPUT_DIR = PROJECT_ROOT / "test_input"
TEST_OUTPUT_DIR = PROJECT_ROOT / "test_outputs" 


# Input for binary matrices
SEEDS_SET = INPUT_DIR / "updated_seedsets_of_interest.pckl"
NON_SEEDS_SET = INPUT_DIR / "updated_non_seedsets_of_interest.pckl"

# Input Files
SEEDS_PICKLE = INPUT_DIR / "seeds_binary_per_patric.pckl"
NON_SEEDS_PICKLE = INPUT_DIR / "non_seeds_binary_per_patric.pckl"


# Testers
# Test input
SEEDS_PICKLE_TEST = TEST_INPUT_DIR / "test_seed_patric.pckl"
NON_SEEDS_PICKLE_TEST = TEST_INPUT_DIR / "test_non_seed_patric.pckl"
COMPACT_METADATA_TEST = TEST_INPUT_DIR / "test_PREGO.csv"

# Test output 
PATRIC_METADATA_JSON_TEST = TEST_OUTPUT_DIR / "patric_ids_metadata.json"
FAILED_PATRIC_IDS_TXT_TEST = TEST_OUTPUT_DIR / "failed_patric_ids.txt"
KEGG_DATA_JSON_TEST = TEST_OUTPUT_DIR / "kegg_flat_entries.json"
FAILED_KEGG_IDS_TXT_TEST = TEST_OUTPUT_DIR / "failed_kegg_ids.txt"
PREGO_ENVIRONMENTS_JSON_TEST = TEST_OUTPUT_DIR / "prego_environments.json"
FAILED_PREGO_IDS_TXT_TEST = TEST_OUTPUT_DIR / "failed_prego_ids.txt"


# Metadata
# Raw Genome Metadata Files
GTDB_METADATA_TSV = METADATA_RAW_DIR / "gtdb_metadata_r207.tsv"
PATRIC_METADATA_JSON = METADATA_RAW_DIR / "patric_ids_metadata.json"
FAILED_PATRIC_IDS_TXT = METADATA_RAW_DIR / "failed_patric_ids.txt"
KEGG_DATA_JSON = METADATA_RAW_DIR / "kegg_flat_entries.json"
FAILED_KEGG_IDS_TXT = METADATA_RAW_DIR / "failed_kegg_ids.txt"
PREGO_ENVIRONMENTS_JSON = METADATA_RAW_DIR / "prego_environments.json"
FAILED_PREGO_IDS_TXT = METADATA_RAW_DIR / "failed_prego_ids.txt"

# Processed Genome Metadata
PREGO_INPUT_CSV = METADATA_DIR / "patric_gtdb_metadata.csv"
COMPACT_METADATA_TSV = METADATA_DIR / "compact_genome_metadata.tsv"
COMPACT_METADATA_WITH_BIOME_TSV = METADATA_DIR / "compact_metadata_with_biome.tsv"

# Compound-level metadata
COMPOUNDS_TSV = METADATA_RAW_DIR / "compounds.tsv"

# Compound-level processed metadata
COMPOUND_SUMMARY_TSV = METADATA_DIR / "compounds_summary.tsv"
BRITE_UNCLASSIFIED_TSV = METADATA_DIR / "brite_unclassified.tsv"
AMINOACIDS_TSV = METADATA_DIR / "common_amino_acid_compounds.tsv"


# Trees
BAC_TREE = GTDB_TREES_DIR/ "bac120_r207.tree"
AR_TREE = GTDB_TREES_DIR/ "ar53_r207.tree"
BAC_TREE_PRUNED = GTDB_TREES_DIR/ "bac120_subset_pruned.tree"
AR_TREE_PRUNED = GTDB_TREES_DIR/ "ar53_subset_pruned.tree"


# Output Files
METABOLIC_POTENTIAL_TSV = OUTPUT_DIR / "metabolic_potential_summary.tsv"
GENOME_PATHWAY_COVERAGE = OUTPUT_DIR / "genome_pathway_coverage.tsv"
