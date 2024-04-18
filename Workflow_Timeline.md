# Introduction

This project addresses the genomic analysis of *A. lyrata* species to understand the genetic diversity, structure, and evolutionary aspects. The analysis includes the identification of population structures, allele frequency estimation, and the investigation of hybrids. The expected outcomes are detailed PCA plots, population structure inferences, and allele frequency distributions across different populations.

# Package Installation

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
  
# Data

## GATK Variant Filtering:
**Language:**
UNIX

### Input File 
- `Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf` : This file contains variant data that will be filtered by the command.

### Output Data:
- `my_filtered_output.vcf`: The VCF file after initial filtering to include specific populations.

### Script used:
NEED TO ADD THIS

## PCA Analysis Using Tomas Script
**Language:**
R

### Input Files
- `my_filtered_output.vcf`
  - This is the same output VCF file used for GATK processing and is reused here.

### Output Files
- The PCA plot is generated and saved as an object within the R environment. 

- Saving the Plot
To save the PCA plot to a file, you can add a `ggsave()` function or similar at the end of the script. This will allow you to save the plot in your desired format and location.

### Script Used
You can view the script used for this analysis from the following link:
- [Tomas Script for PCA](https://github.com/thamala/polySV/blob/main/est_cov_pca.r)


## PCA Analysis Using Adegenet
**Language:**
R

### Input Files
- `my_filtered_output.vcf`
  - This is the same output VCF file used for GATK processing and is reused here.

### Output Files
- PCA plots saved as PDF: `PCA_all_SNPs_ax12_1K_less.pdf` 

### Script Used
ADD HERE
  - Provided by Dr. Levi Yant


### PCA Method Comparison
Both scripts process the same VCF file but employ different methods tailored to specific aspects of the data:
- **Tomas Script:** Focuses on filtering SNPs based on allele frequency, aiming for a refined dataset for PCA.
- **Adegenet Script:** Optimized for speed and efficiency in PCA calculations, particularly effective for large-scale genomic data analysis.


## Detailed PCA and NJ Tree Analysis Using Adegenet
**Language:**
R

### Input Files
- `my_filtered_output.vcf`
  - This is the same output VCF file used for GATK processing and is reused here.

### Output Files
- **Graphical Outputs:**
  - PCA plots saved as a PDF: `PCA_populations_location.pdf`
  - Ploidy-differentiated PCA plots saved as a PDF: `Ana_v_PCA_all_ploidycol_SNPs_ax12_1K_less.pdf`
- **Distance Matrices:**
  - Individual distance matrix in Phylip format: `aa.indiv_Neis_distance_4ds.phy.dst`
  - Population distance matrix in Phylip format: `aa.pops_Neis_distance_4ds.phy.dst`
- **Tree Files:**
  - NJ tree of populations: `NJ.tree_populations_4ds.tre`
  - NJ tree based on distance analysis: `NJ.distance_tree_outgroups.tre`  

### Script Used
ADD HERE
  - Provided by Dr. Levi Yant

### Additional Analysis Info
This script extends the functionalities seen in previous PCA scripts by incorporating geographic distance calculations and Nei's genetic distance. It facilitates the creation of Neighbour Joining (NJ) trees and distance-based analyses, offering a deeper exploration of relationships and distances between samples or populations. This can reveal more complex patterns than those discerned by traditional phylogenetic trees.

The `.dst` file can then be used for SplitsTree, to explore relationships and distances between your samples or populations in a way that might reveal more complex patterns than traditional phylogenetic trees.

### Tool for Further Analysis
- **SplitsTree:** Use the following link to download and utilize SplitsTree for detailed analysis of the distance matrices generated:
  - [Download SplitsTree](https://uni-tuebingen.de/en/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/algorithms-in-bioinformatics/software/splitstree/)


# Outlier Removal Process

## GATK Variant Filtering
**Language:**
UNIX

### Input File 
- `Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf` : This file contains variant data that will be filtered by the command.

### Output Data:
- `my_filtered_output_no_outliers.vcf`: The VCF file after initial filtering to include specific populations.

### Script Used
ADD HERE



## PCA Analysis Using Tomas Script Post-Outlier Removal
**Language:**
R

### Input Files
- `my_filtered_output_no_outliers.vcf`
  - This is the same output VCF file used for GATK processing and is reused here.

### Output Files
- The PCA plot is generated and saved as an object within the R environment. 

- Saving the Plot
To save the PCA plot to a file, you can add a `ggsave()` function or similar at the end of the script. This will allow you to save the plot in your desired format and location.

### Script Used
You can view the script used for this analysis from the following link:
- [Tomas Script for PCA](https://github.com/thamala/polySV/blob/main/est_cov_pca.r)



## PCA Analysis Using Adegenet Post-Outlier Removal
**Language:**
R

### Input Files
- `my_filtered_output_no_outliers.vcf`
  - This is the same output VCF file used for GATK processing and is reused here.

### Output Files
- PCA plots saved as PDF: `PCA_all_SNPs_ax12_1K_less_2.pdf` (changed name from script to be able to differentiate)

### Script Used
ADD HERE
  - Provided by Dr. Levi Yant


## Detailed PCA and NJ Tree Analysis Using Adegenet Post-Outlier Removal
**Language:**
R

### Input Files
- `my_filtered_output_no_outliers_2.vcf`
  - This is the same output VCF file used for GATK processing and is reused here.

### Output Files
changed default names from script by adding '_2'
- **Graphical Outputs:**
  - PCA plots saved as a PDF: `PCA_populations_location_2.pdf`
  - Ploidy-differentiated PCA plots saved as a PDF: `Ana_v_PCA_all_ploidycol_SNPs_ax12_1K_less_2.pdf`
- **Distance Matrices:**
  - Individual distance matrix in Phylip format: `aa.indiv_Neis_distance_4ds.phy_2.dst`
  - Population distance matrix in Phylip format: `aa.pops_Neis_distance_4ds.phy_2.dst`
- **Tree Files:**
  - NJ tree of populations: `NJ.tree_populations_4ds_2.tre`
  - NJ tree based on distance analysis: `NJ.distance_tree_outgroups_2.tre`

 ### Additional Information
- The `.dst` files are particularly used for SplitsTree analysis, facilitating the exploration of more complex genetic relationships.


### Script Used
ADD HERE
  - Provided by Dr. Levi Yant


## Allele Frequency Analysis
**Language:** UNIX

### Input Files 
#### With Outliers
- `my_filtered_output.vcf`
  - The VCF file containing variant calls including outliers.
- `with_individuals_population_ids.txt`
  - A text file mapping individuals to their populations including outliers.

#### Without Outliers
- `my_filtered_output_no_outliers.vcf`
  - The VCF file after outliers have been removed.
- `without_individuals_population_ids_updated.txt`
    - Updated population file reflecting the removal of outliers.

### Output Files 
- Allele Frequency Output (With Outliers): `freq_with`
  - Contains allele frequencies calculated from the dataset that includes outliers.
- Allele Frequency Output (Without Outliers): `freq_without`
  - Contains allele frequencies calculated from the dataset after removing outliers.

### Scripts Used
- **Script for Allele Frequencies:** `poly_freq.c`
  - A custom C program that calculates allele frequencies from a VCF file.
- **Compiler:** gcc
  - The GCC compiler is used to compile the C program.
- **Script for Population Files:** Script to make tab delimited text file (with and without)
  - This script prepares population data files for frequency calculation.

### Additional Information
The process for calculating allele frequencies uses the custom C program `poly_freq.c`, which operates with command-line parameters. This program requires the installation of gcc for compilation. The process provides allele frequency outputs for datasets with and without outliers, allowing for comparative analysis of genetic variation influenced by population structure.



## Allele Frequency Histogram Analysis
**Language:** R

### Input Files 

  - `freq_with.txt`: This text file contains allele frequency data for the population including outliers. The data include columns for CHROM (chromosome), POS (position), and AF (allele frequency).
  - `freq_without.txt`: Similar to `freq_with.txt`, this file contains allele frequency data for the population excluding outliers, with the same columns for chromosome, position, and allele frequency.

### Output Files 
- **Histogram of Hybrid Simulation:**
  - `Histogram_of_Hybrid_Simulation.pdf`
  - This PDF document contains a histogram representing the mean allele frequencies across the combined datasets. The histogram provides a visual comparison of genetic variation with and without the influence of outliers.
  - Note: This file is created and saved in the current working directory. You can adjust the file name and path as needed using the `pdf()` function call in the R script.

### Script Used














