from pathlib import Path


# Project Root
PROJECT_ROOT = Path(__file__).resolve().parent


# Directory Structure
INPUT_DIR = PROJECT_ROOT / "input"
OUTPUT_DIR = PROJECT_ROOT / "output"
TEST_DIR = PROJECT_ROOT / "test_input" 
PLOTS_DIR = OUTPUT_DIR / "plots"
MAPS_DIR = PROJECT_ROOT / "maps"

METADATA_DIR = PROJECT_ROOT / "metadata"
METADATA_RAW_DIR = METADATA_DIR / "raw"

SCRIPTS_DIR = PROJECT_ROOT / "scripts"
UTILS_DIR = PROJECT_ROOT / "utils"

GTDB_TREES_DIR = PROJECT_ROOT / "gtdb_trees"


# Input Files
SEEDS_PICKLE = INPUT_DIR / "seeds_binary_per_patric.pckl"
#SEEDS_PICKLE = TEST_DIR / "test_patric.pckl"
NON_SEEDS_PICKLE = INPUT_DIR / "non_seeds_binary_per_patric.pckl"
SEEDS_SET = INPUT_DIR / "updated_seedsets_of_interest.pckl"
NON_SEEDS_SET = INPUT_DIR / "updated_non_seedsets_of_interest.pckl"


# Metadata Files
GTDB_METADATA = METADATA_RAW_DIR / "gtdb_metadata_r207.tsv"
PATRIC_METADATA_JSON = METADATA_RAW_DIR / "patric_ids_metadata.json"
FAILED_PATRIC_IDS_TXT = METADATA_RAW_DIR / "failed_patric_ids.txt"
COMPACT_METADATA = METADATA_DIR / "compact_genome_metadata.csv"
COMPACT_METADATA_WITH_BIOME = METADATA_DIR / "compact_metadata_with_biome.tsv"

REACTIONS_OF_INTEREST = METADATA_DIR / "reactions_of_interest.tsv"


# Compound-level metadata
COMPOUNDS_TSV = METADATA_RAW_DIR / "compounds.tsv"
REACTIONS_TSV = METADATA_RAW_DIR / "reactions.tsv"
PATHWAYS_TSV = METADATA_RAW_DIR / "KEGG.pathways"
KEGG_DATA = METADATA_RAW_DIR/"kegg_flat_entries.json"
COMPOUND_SUMMARY_TSV = METADATA_DIR / "compounds_summary.tsv"
CPDs_KEGG_DATASET_TSV = METADATA_DIR / "cpds_KEGG_dataset.tsv"
#CPDs_KEGG_DATASET_TSV = TEST_DIR / "test_KEGG.tsv"


# Maps 
MODULE_REACTION_MAPPING = MAPS_DIR / "module_reaction.tsv"
CPD_MODULE_MAPPING = MAPS_DIR/ "KEGG_cpd_module.tsv"
CPD_PATHWAY_MAPPING = MAPS_DIR/ "cpd_map.tsv"
CPD_REACTION_MAPPING = MAPS_DIR / "cpd_reaction.tsv"
PATHWAY_REACTION_MAPPING = MAPS_DIR / "pathway_reaction.tsv"


# Trees
BAC_TREE = GTDB_TREES_DIR/ "bac120_r207.tree"
AR_TREE = GTDB_TREES_DIR/ "ar53_r207.tree"
BAC_TREE_PRUNED = GTDB_TREES_DIR/ "bac120_subset_pruned.tree"
AR_TREE_PRUNED = GTDB_TREES_DIR/ "ar53_subset_pruned.tree"


# Output Files
METABOLIC_POTENTIAL_1 = OUTPUT_DIR / "metabolic_potential_summary_per_mbp.csv"
#METABOLIC_POTENTIAL_1 = TEST_DIR / "test_PREGO.csv"
GENOME_PATHWAY_COVERAGE = OUTPUT_DIR / "genome_pathway_coverage.tsv"
AMINOACIDS = OUTPUT_DIR / "common_amino_acid_compounds.tsv"
UNIQUE_REACTIONS = OUTPUT_DIR / "unique_kegg_reactions.tsv"

PREGO_ENVIRONMENTS_JSON = METADATA_RAW_DIR / "prego_environments.json"
FAILED_PREGO_IDS_TXT = METADATA_RAW_DIR / "failed_prego_ids.txt"