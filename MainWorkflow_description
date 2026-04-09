## Main workflow

### 1. Generate input interaction/domain files and run Chrom3D

The main pipeline is provided in:

- `NCHG_Chrom3D_modeling.txt`

This script performs the following steps:

1. Converts `.hic` files to `.mcool`
2. Computes expected cis/trans contacts using **NCHG**
3. Processes TAD/domain BED files
4. Builds domain cartesian products for NCHG testing
5. Computes significant cis/trans interactions
6. Filters interactions by **FDR** and **log-ratio**
7. Defines peripheral/LAD-like regions from low compartment-score bins
8. Builds `.gtrack` files for Chrom3D
9. Adds LAD/periphery annotation
10. Runs multiple Chrom3D simulations in parallel

#### Key outputs
- `cool_files/*.mcool`
- `h5/*`
- `parquet/*`
- `sig_interactions/*`
- `gtrack/*`
- `sim_all_lr2_2mil/*.cmm`
- `sim_all_lr2_2mil/*.log`

---

### 2. Build Chrom3D genome graph inputs

#### `makeGtrack.py`
Creates a Chrom3D `.gtrack` file from:

- a significant interaction BED file
- a genome/domain BED file

This script converts significant pairwise interactions into Chrom3D graph edges.

#### `lad_gtrack.py`
Adds periphery/LAD information to an existing `.gtrack` file by overlapping domains with a BED-like set of peripheral regions.

This is used to mark bins/beads as nuclear periphery-associated before Chrom3D simulation.

---

### 3. Monitor Chrom3D optimization

#### `plot_loss_score.R`
Plots Chrom3D loss trajectories across multiple simulation replicates for a given group of samples.

Useful for:
- checking optimization stability
- comparing convergence across age groups or samples
- selecting well-behaved simulation sets

---

### 4. Assign Chrom3D beads to subcompartments

#### `model_to_subcom.py`
Maps Chrom3D beads to subcompartment labels by overlapping model coordinates with a BED file of subcompartments.

It collapses subcompartment labels into a simplified 8-state scheme:

- `A3`, `A2`, `A1`, `A0`
- `B0`, `B1`, `B2`, `B3`

#### Output
- bead-level subcompartment annotation table

---

### 5. Detect spatial SNP hotspot regions

#### `detect_most_freq_SNP_domains.py`
Identifies Chrom3D beads/domains overlapping SNP coordinates, counts trait enrichments, and selects top SNP-dense beads.

The script:
- overlaps model beads with SNP coordinates
- counts **immune-related** and **neurodegenerative** SNP burden
- keeps top-ranked beads
- computes pairwise Euclidean distances among selected beads

#### Outputs
- `*_top30.tsv` for bead-to-bead distance matrices
- `*_top30_beads.tsv` for selected bead annotations

---

### 6. Cluster SNP-enriched beads in 3D space

#### `SNP_heatmaps.R`
Uses pairwise Euclidean distance matrices from top SNP-enriched beads to:

- cluster beads hierarchically
- choose an approximate optimal cluster number using silhouette width
- generate clustered heatmaps
- export bead-to-cluster assignments

#### Outputs
- cluster heatmaps (`.pdf`)
- cluster assignments (`.tsv`)

---

### 7. Summarize cluster composition and radial positioning

#### `snp_radial_dist_per_cluster.R`
Computes distance of SNP-enriched beads from the nuclear center and compares radial positioning across clusters.

Useful for testing whether specific SNP-defined spatial clusters are:

- more central
- more peripheral
- structurally segregated

#### `cluster_gwas_share.R`
Calculates normalized GWAS trait composition per cluster and visualizes cluster-level trait proportions.

Example traits in the current script include:

- `AD`
- `PD`
- `MS`
- `IBD`

Normalization is performed relative to trait set size.

#### `complex_upset_plots.R`
Creates UpSet-style summaries of trait overlap across cell types for clustered SNP-associated regions.

Useful for identifying:

- shared versus cell type-specific disease-associated spatial hotspots
- trait intersection structure across samples

---

### 8. Functional enrichment of spatial clusters

#### `GO_per_cluster.R`
Performs GO Biological Process enrichment for genes overlapping clustered regions.

The script:

- maps cluster intervals to genes using **biomaRt**
- groups genes by cluster
- runs `compareCluster()` from **clusterProfiler**
- outputs dotplots and Excel tables

#### Outputs
- `GO_<sample>.pdf`
- `GO_<sample>.xlsx`

#### `GO_per_cluster_summary_plot.R`
Aggregates GO enrichment results across samples and visualizes shared cluster-level themes using scatterpie summaries.

This is useful for highlighting whether spatial clusters are predominantly associated with:

- immune processes
- neurodegenerative processes
- mixed functional programs

---

### 9. Subcompartment composition of spatial clusters

#### `subcom_share_per_cluster.R`
Merges cluster assignments with bead-level subcompartment labels and plots the subcompartment composition of each cluster.

Useful for asking whether specific spatial SNP clusters are enriched for:

- active compartment-like regions
- inactive/repressive compartment-like regions
- peripheral B3-like states

#### `sub_compartment_radial_positioning_analysis.R`
Summarizes radial positioning of subcompartments across cell types and age groups.

This script computes median radial position for each subcompartment across repeated Chrom3D simulations and visualizes:

- subcompartment organization by **cell type**
- subcompartment organization by **age group**

---

### 10. Network topology of interaction structures

#### `networks_properties.R`
Builds graph representations of significant domain-domain interactions and summarizes global network properties, including:

- community number
- global clustering
- average local clustering
- cis share
- trans share

Also exports a summary table of graph properties.

#### `networks_properties_complete_table.R`
A more detailed graph-level summary script that additionally computes node-level or graph-level metrics such as:

- degree
- closeness
- betweenness
- eigenvector centrality
- maximum clique size
- average path length
- diameter
- assortativity

Use this version when a more complete topology table is needed for downstream statistical analysis or supplementary materials.

---

## Suggested analysis order

For a standard end-to-end run, the recommended order is:

1. **Run NCHG and prepare Chrom3D inputs**
   - `NCHG_Chrom3D_modeling.txt`
   - `makeGtrack.py`
   - `lad_gtrack.py`

2. **Run Chrom3D simulations**
   - Chrom3D output: `.cmm`, `.log`

3. **Inspect optimization**
   - `plot_loss_score.R`

4. **Run whole-network analyses**
   - `networks_properties.R`
   - `networks_properties_complete_table.R`

5. **Summarize genome organization**
   - `sub_compartment_radial_positioning_analysis.R`

6. **Annotate beads**
   - `model_to_subcom.py`
   - `detect_most_freq_SNP_domains.py`

7. **Cluster SNP-enriched beads**
   - `SNP_heatmaps.R`

8. **Perform downstream cluster analyses**
   - `snp_radial_dist_per_cluster.R`
   - `cluster_gwas_share.R`
   - `complex_upset_plots.R`
   - `subcom_share_per_cluster.R`
   - `GO_per_cluster.R`
   - `GO_per_cluster_summary_plot.R`

---

## Software requirements

### Core tools
- [Chrom3D](https://github.com/Chrom3D/Chrom3D)
- [NCHG](https://github.com/paulsengroup/NCHG) or local compiled installation
- `hictk`
- `bedtools`
- GNU `parallel`

### Python packages
- `pandas`
- `numpy`
- `bioframe`
- `scipy`

### R packages
- `tidyverse`
- `dplyr`
- `tidyr`
- `ggplot2`
- `pheatmap`
- `cluster`
- `RColorBrewer`
- `ComplexUpset`
- `scatterpie`
- `biomaRt`
- `clusterProfiler`
- `org.Hs.eg.db`
- `openxlsx`
- `igraph`
- `ggpubr`

---

## Input data expected

Depending on the script, the repository expects some combination of the following:

- `.hic` or `.mcool` Hi-C files
- TAD/domain BED files
- compartment/subcompartment BED or BEDGRAPH files
- SNP coordinate tables
- Chrom3D `.cmm.tsv` model files
- Chrom3D `.log` files
- cluster assignment tables
- bead-level coordinate tables

Many scripts currently use placeholder paths such as `input_dir <- ""` or `input_dir = "dir"`, so these should be edited before running.

---

