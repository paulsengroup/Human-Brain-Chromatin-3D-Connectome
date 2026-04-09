# 3D Genome Modeling and Downstream Analysis

A collaborative repository for **Chrom3D-based 3D genome modeling** and downstream spatial analyses integrating:

- significant chromatin interactions from **NCHG**
- **TAD/domain-based** genome graph construction
- **LAD/periphery-aware** Chrom3D simulations
- **subcompartment assignment**
- **GWAS/SNP spatial hotspot** analyses
- **cluster-level functional annotation**
- **network topology** summaries across cell types and age groups

This repository is designed for user-friendly, modular analysis of **cell type-resolved 3D genome architecture** across biological conditions.

---

## Overview

This workflow starts from processed Hi-C data and domain annotations, builds Chrom3D-compatible genome graphs, runs multiple 3D simulations, and performs downstream analyses on:

- 3D spatial proximity of disease-associated loci
- cluster structure of top SNP-enriched beads
- radial positioning and nuclear organization
- subcompartment composition
- GO enrichment of spatial clusters
- graph/network properties of interaction structures

The codebase combines:

- **bash** for end-to-end preprocessing and modeling
- **Python** for interval overlap, annotation, and formatting
- **R** for clustering, enrichment, summaries, and visualization

---

## Repository contents

```text
.
├── NCHG_Chrom3D_modeling.txt
├── makeGtrack.py
├── lad_gtrack.py
├── model_to_subcom.py
├── detect_most_freq_SNP_domains.py
├── SNP_heatmaps.R
├── snp_radial_dist_per_cluster.R
├── cluster_gwas_share.R
├── complex_upset_plots.R
├── GO_per_cluster.R
├── GO_per_cluster_summary_plot.R
├── subcom_share_per_cluster.R
├── sub_compartment_radial_positioning_analysis.R
├── networks_properties.R
├── networks_properties_complete_table.R
└── plot_loss_score.R
