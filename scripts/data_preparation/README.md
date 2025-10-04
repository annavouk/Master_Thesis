# Data Preparation and Exploratory Overview

This folder contains scripts for preparing genome- and compound-level metadata 
for downstream analyses. It integrates raw resources from GTDB, PATRIC, PREGO, 
ModelSEED, and KEGG into unified, analysis-ready tables.

Input pickles were downloaded from [microbetag](https://zenodo.org/records/17172592)'s Zenodo repository (seeds_per_genome.pkl.gz and nonseeds_per_genome.pkl.gz). 

GTDB metadata [r207](https://data.gtdb.ecogenomic.org/releases/release207/207.0/) (ar53_metadata_r207.tar.gz and bac120_metadata_r207.tar.gz) was downloaded and merged in gtdb_metadata_r207.tsv.

```console
    head -n 1 metadata/raw/ar53_metadata_r207.tsv > metadata/raw/gtdb_metadata_r207.tsv
    tail -n +2 -q metadata/raw/ar53_metadata_r207.tsv metadata/raw/bac120_metadata_r207.tsv >> metadata/raw/gtdb_metadata_r207.tsv
```

compounds.tsv was downloaded with wget from https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/refs/heads/master/Biochemistry/compounds.tsv


## Data Preparation Scripts

- **fetch_metadata.py**  
  Generic metadata fetcher for KEGG compounds, PATRIC genomes, and PREGO 
  environmental annotations.

- **merge_patric_gtdb_metadata.py**  
  Merges genome metadata from GTDB (release r207) and PATRIC API outputs into a 
  compact summary table. Adds environmental annotations (PREGO) and genome 
  quality metrics (CheckM).

- **merge_cpds_reactions_modules_pathways_ontology.py**  
  Integrates ModelSEED compounds with KEGG identifiers, reactions, modules, 
  pathways, and BRITE ontology (from bulk-downloaded KEGG flat files).

- **maps_envs_to_biomes.py**  
  Groups PREGO environment annotations into broad biome categories using 
  controlled regex rules.

### Inputs

- input/seeds_per_genome.pkl 
- input/nonseeds_per_genome.pkl 
- metadata/raw/gtdb_metadata_r207.tsv (GTDB genomes)  
- metadata/raw/patric_metadata.json (PATRIC genomes)  
- metadata/raw/prego_envs.json (PREGO annotations)  
- metadata/raw/modelseed_biochemistry.tsv (ModelSEED compounds)  
- metadata/raw/kegg_flat_entries.json (KEGG compound entries)

### Outputs

- metadata/compact_metadata.tsv (intermediate file)
- metadata/compact_metadata_with_biome.tsv
  Unified genome-level metadata table with taxonomy, genome quality, and biome annotations.

- metadata/compound_summary.tsv
  Unified compound-level metadata with KEGG IDs, reactions, modules, pathways, and BRITE ontology.


## Exploratory Overview Scripts

- **overview_genome_taxonomy.py**  
  Summarizes the taxonomic distribution of the genome dataset and produces plots:  
  - Barplot of unique taxa per rank  
  - Pie chart of genomes per domain  
  - Barplots of top 10 phyla and genera  

- **overview_compounds.py**  
  Summarizes the distribution of metabolic compounds and produces plots:  
  - Venn diagram of seed vs non-seed overlap  
  - Histograms of KEGG modules and pathways per compound  
  - Barplots of top KEGG pathways and modules  
  - Barplot of BRITE ontology distribution  

### Inputs

- input/seeds_per_genome.pkl  
- input/nonseeds_per_genome.pkl 
- metadata/compact_metadata_with_biome.tsv (genome-level metadata)  
- metadata/compound_summary.tsv (compound-level metadata)  

### Outputs

Plots are exported to outputs/plots/exploratory_data_overview/


## Notes

- Scripts here provide descriptive summaries only.  
- Results serve as a foundation for downstream comparative and statistical analyses.  
