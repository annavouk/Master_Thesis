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
AMINOACIDS_PLOTS_DIR = PLOTS_DIR / "aminoacids_analysis"
INTERACTION_ANALYSIS_PLOTS_DIR = PLOTS_DIR / "interaction_analysis"

METADATA_DIR = PROJECT_ROOT / "metadata"
METADATA_RAW_DIR = METADATA_DIR / "raw"

SCRIPTS_DIR = PROJECT_ROOT / "scripts"

UTILS_DIR = PROJECT_ROOT / "utils"

TEST_INPUT_DIR = PROJECT_ROOT / "test_input"
TEST_OUTPUT_DIR = PROJECT_ROOT / "test_outputs" 


# Input Files
SEEDS_PICKLE = INPUT_DIR / "seeds_per_genome.pkl"
NON_SEEDS_PICKLE = INPUT_DIR / "nonseeds_per_genome.pkl"


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
KEGG_PATHWAYS_TSV = METADATA_RAW_DIR / "KEGG.pathways"

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


# Output Files
METABOLIC_POTENTIAL_TSV = OUTPUT_DIR / "metabolic_potential_summary.tsv"

KRUSKAL_BIOMES_SUBSAMPLES_TSV = OUTPUT_DIR / "kruskal_subsampling.tsv"
KRUSKAL_BIOMES_SUMMARY_TSV = OUTPUT_DIR / "kruskal_summary.tsv"

DUNN_BIOMES_SUBSAMPLES_TSV = OUTPUT_DIR / "dunn_posthoc_subsampling.tsv"
DUNN_BIOMES_SUMMARY_TSV = OUTPUT_DIR / "dunn_posthoc_summary.tsv"

OUTLIERS_L2H_EDGES_CSV = OUTPUT_DIR / "ratio_extreme_low2high_edgelist.csv"
OUTLIERS_H2L_EDGES_CSV = OUTPUT_DIR / "ratio_extreme_high2low_edgelist.csv"
CYANO_HALO_EDGES_CSV = OUTPUT_DIR / "cyano_halo_all_edgelist.csv"
