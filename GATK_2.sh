#Refine Variant Selection
## This command filters the VCF file to include only specified samples while excluding outliers identified in previous analyses. Adjust the paths to your local directories and ensure the correct sample names are used.
## The `-sn` flag includes specific samples in the variant filtering process. Only samples listed will be included in the output VCF, ensuring that data from unwanted samples or identified outliers are excluded. Make sure to review and update the sample names to reflect those in your study that meet the current criteria for inclusion.
## The `-O` option specifies the output location for the newly created VCF file. Change the output file path as needed to match your directory structure or naming conventions. This file will contain only the filtered variants from the selected samples.


### Activate GATK Environment
conda activate gatk4


gatk SelectVariants \
-V /Users/yaz/Desktop/r3/Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf \
-sn KAG-01tl -sn KAG-02tl -sn KAG-03tl -sn KAG-04tl \
-sn LIC-01tl -sn LIC-02tl -sn LIC-03tl \
-sn MOD-01tl -sn MOD-02tl -sn MOD-03tl -sn MOD-04tl \
-sn FRE-01tl -sn FRE-02tl -sn FRE-03tl -sn FRE-04tl -sn FRE-05tl -sn FRE-07tl -sn FRE-06tl -sn FRE-08tl \
-sn OCH-01tl -sn OCH-02tl -sn OCH-03tl -sn OCH-04tl -sn OCH-06tl -sn OCH-07tl -sn OCH-08tl -sn OCH-05tl \
-sn HAB-01tl -sn HAB-02tl -sn HAB-03tl \
-sn PIZ-03dl -sn PIZ-04dl -sn PIZ-05dl -sn PIZ-06dl -sn PIZ-08dl -sn PIZ-09dl -sn PIZ-11dl \
-sn VLH-02dl -sn VLH-03dl -sn VLH-05dl -sn VLH-06dl -sn VLH-07dl -sn VLH-01dl \
-sn KEH-01tl -sn KEH-04tl -sn KEH-05tl -sn KEH-10tl -sn KEH-09tl -sn KEH-07tl -sn KEH-06tl -sn KEH-08tl \
-sn BZD-01tl -sn BZD-02tl -sn BZD-03tl -sn BZD-04tl -sn BZD-05tl -sn BZD-06tl -sn BZD-07tl -sn BZD-08tl \
-O /Users/yaz/Desktop/r3/my_filtered_output_no_outliers.vcf


### Dectivate GATK Environment
conda deactivate



