# Core Analysis

This directory contains the main scripts implementing the thesis analyses of microbial metabolic potential and interactions.
Analyses include genome-level metabolic potential metrics, biome-level comparisons, seed node frequencies, amino acid auxotrophy/prototrophy patterns, and provider–receiver interactions.

The directory is organized into two subfolders:
- **analysis/** – core computations (metrics, biome-level tests, auxotrophy, provider-receiver)  
- **visualization/** – generation of plots and figures corresponding to these analyses


## Analysis Scripts

- **metabolic_potential.py**
Calculates total Seeds, total non-Seeds, their ratio per genome.  

- **make_kruskal_wallis_summary_subsampling.py** 
Computes summary statistics per Metric for Kruskal-Wallis results of 100 subsamples of original dataset. 

- **make_dunn_summary_subsampling.py**
Computes summary statistics per Metric for dunn's post-hoc results of 100 subsamples of original dataset.

- **provider_receiver.py** 
Builds feeding edges (provider-receiver interactions)

### Inputs

- input/seeds_per_genome.pkl  
- input/nonseeds_per_genome.pkl 
- metadata/compact_genome_metadata.tsv
- output/kruskal_subsampling.tsv
- output/dunn_posthoc_subsampling.tsv"
- metadata/compact_metadata_with_biome.tsv
- output/metabolic_potential_summary.tsv

### Outputs

Output files are saved to outputs/.



## Visualization Scripts

- **metabolic_potential_visualization.py**  
Performs exploratory and statistical analysis on genome-wide metabolic potential
(measured as the ratio of seed to non-seed compounds per genome), and visualizes the data using histograms, scatterplots, violin plots and barplots and assess differences between taxonomic groups with Kruskal-Wallis and Dunn's post-hoc tests.

- **metabolic_potential_by_biome.py** 
Compares genome-level metabolic potential across soil, marine and freshwater biomes. Generates boxplots and stacked barplot.

- **seed_node_frequency.py**
Analyzes the frequency of seed nodes (essential metabolites) across genomes and generates histogram and barplot.

- **cpd_aminoacids.py**
Examines amino acid Auxotrophy & Prototrophy profiles across genomes. Generates barplot, histograms, and stacked barplots.

- **provider_receiver_visualization.py** 
Processes edge lists of microbial interactions, computes top providers/receivers, and plots heatmaps. 

### Inputs

- output/metabolic_potential_summary.tsv
- metadata/compact_metadata_with_biome.tsv
- input/seeds_per_genome.pkl  
- input/nonseeds_per_genome.pkl 
- metadata/compounds_summary.tsv
- metadata/compact_genome_metadata.tsv
- metadata/common_amino_acid_compounds.tsv
- output/ratio_extreme_low2high_edgelist.csv
- output/ratio_extreme_high2low_edgelist.csv
- output/cyano_halo_all_edgelist.csv

### Outputs

Plots are exported to:
- outputs/plots/metabolic_potential
- outputs/plots/biomes_analysis
- outputs/plots/seed_node_frequency
- outputs/plots/aminoacids_analysis
- outputs/plots/interaction_analysis
