# Introduction

This project addresses the genomic analysis of *A. lyrata* species to understand the genetic diversity, structure, and evolutionary aspects. The analysis includes the identification of population structures, allele frequency estimation, and the investigation of hybrids. The expected outcomes are detailed PCA plots, population structure inferences, and allele frequency distributions across different populations.

## Package Installation

Install specific versions of the packages used in this analysis:

**Install the following R packages:**
- `vcfR=1.15.0`
- `ggplot2=3.4.4`
- `ggrepel=0.9.5`
- `adegenet=2.1.10`
- `pegas=1.3`
- `StAMPP=1.6.3`
- `ade4=1.7.22`
- `MASS=7.3.60`
- `adegraphics=1.0.21`
- `geosphere=1.5.18`
- `dplyr=1.1.4`
- `tidyverse=2.0.0`
- `ape=5.7.1`

**Other Tools:**
- `gatk4= 4.5.0.0`
- `gcc 14.0.0` (install using Homebrew)
  
## Script Files

- `est_cov_pca.r`: An R script file for estimating the covariance matrix and performing PCA.
- `Easy_adegenet_script_and_modified_functions.R`: An R script file for advanced PCA plotting and generating `.dst` files for further analysis in SplitsTree.
- `poly_freq.c`: C source code for calculating allele frequencies from a VCF file.
  
## Data

### GATK:
language: UNIX

**Input File** 
- VCF File: `my_filtered_output.vcf` - The VCF file after initial filtering to include specific populations.

**Output Data:**
- `my_filtered_output.vcf`: The VCF file after initial filtering to include specific populations.

### Script used:
NEED TO ADD THIS

## PCA Analysis Using Tomas Script

### Input Files
- **VCF File:** `my_filtered_output.vcf`
  - This is the same output VCF file used for GATK processing and is reused here for outlier removal.

**Output Data:**
The PCA plot is generated and saved as an object within the R environment. 

-Saving the Plot
To save the PCA plot to a file, you can add a `ggsave()` function or similar at the end of the script. This will allow you to save the plot in your desired format and location.

### Script Used
You can view and download the script used for this analysis from the following link:
- [Tomas Script for PCA](https://github.com/thamala/polySV/blob/main/est_cov_pca.r)














