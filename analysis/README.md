# Deep Mutational Screening Analysis

This repository contains a series of R Markdown notebooks detailing the analysis of deep mutational screening data for Plasmodium KRS. The goal of this project is to identify potential resistance mutants against wildtype inhibitors through a comprehensive analysis including initial data processing, validation, normalization, and differential expression analysis using edgeR.

## Notebooks Overview

### 01_mutlib-initial-processing

This notebook focuses on the initial processing of the mutation library data. It includes steps for reading datasets, parsing mutation names, translating codons to amino acids, and preparing the data for further analysis.

- [View Notebook](01_mutlib-initial-processing.md)

### 02_mutlib-validation

The second notebook validates the processed data, comparing pre- and post-selection libraries, analysing codon bias, and performing basic statistical analyses to ensure data integrity.

- [View Notebook](02_mutlib-validation.md)

### 03_mutlib-normalisation

In this notebook, we perform normalisation of the mutation library data, adjusting counts by total library size and relative to the pre-selection library, and identify putative resistance alleles.

- [View Notebook](03_mutlib-normalisation.md)

### 04_mutlib-edgeR

The final notebook uses edgeR to calculate post-treatment enrichment, treating each drug as a replicate for significance calculation, and includes visualisation of the results through volcano plots and mapping to 3D structures.

- [View Notebook](04_mutlib-edgeR.md)

## Results and Figures

The analysis generates a series of results and figures, which are stored in the `results/` and `figures/` directories, respectively. These include tables summarizing the enrichment of mutations post-treatment, volcano plots highlighting significant mutations, and interactive plots for detailed examination of specific sites.

- **Results Directory**: Contains CSV files with detailed analysis results, including treatment enrichments and summary tables.
- **Figures Directory**: Includes static and interactive visualizations of the analysis, such as density plots, volcano plots, and ridge plots, alongside mappings to 3D protein structures.

## Getting Started

To run these notebooks, you will need R and several packages including `dplyr`, `tidyr`, `ggplot2`, `edgeR`, and others. Installation instructions for R and these packages can be found at [CRAN](https://cran.r-project.org) and [Bioconductor](https://bioconductor.org).

### Data

The data used in this analysis is available in the root-level `data/` directory. This includes the pre- and post-selection mutation libraries, as well as the wildtype sequence and the library design.

