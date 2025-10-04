# Microbial Metabolic Interaction Potential

This repository contains the code and data processing pipelines developed for the MSc thesis:  
**"Investigation of microbial metabolic interaction potential across the GTDB representative genomes"**.

The project systematically analyzes metabolic potential and interaction potential across **33,755 GTDB representative genomes** by integrating **PATRIC, GTDB, PREGO, ModelSEED, and KEGG** resources.  
The aim is to quantify biosynthetic autonomy versus dependence, detect auxotrophy patterns, and investigate taxonomic and biome-level structuring of metabolic complementarities.


## Repository Structure
```console
.
├── README.md                        # Repository overview and usage
├── config.py                        # Centralized configuration of paths
├── gtdb_trees                       # GTDB phylogenetic trees (bac120, ar53)
├── input                            # Binary matrices (Seeds, Non-Seeds)
├── metadata                         # Raw and processed metadata (taxonomy, genome size, biomes)
├── notes.md                         # Development notes
├── output                           # Analysis results and figures
├── requirements.txt                 # Dependency versions for reproducibility
├── scripts                          # Data preparation, analysis and visualization scripts
├── test_input                       # Test data
├── test_outputs                     # Test results
└── utils                            # Custom helper functions
```

## Warning

**Metadata preparation**

 The metadata files are compressed (.gz). Before running any scripts, please decompress them (e.g., using gunzip) and place the resulting .tsv files in the metadata/ folder.

For example:

```console
gunzip metadata/compact_metadata.tsv.gz
# or keep the original .gz file
gunzip -k metadata/compact_metadata.tsv.gz
```
