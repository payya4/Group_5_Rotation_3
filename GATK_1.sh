
# Ensure that you are using a compatible version of GATK for this script. This script was tested with GATK version 4.5.0.0
# MAC OS version 12.6.5
# Terminal version 2.12.7


# Environment activation                                                                                                                         

conda activate gatk4                                                    

# Ensure that the `gatk4` environment has been properly installed and configured before attempting to activate it. If not already set up, refer to the GATK installation guide.
# Replace the input VCF file path and the output file path with the paths relevant to your project. The sample names (`-sn`) should also be customized to match the samples you wish to include in the analysis.
# The `-sn` options specify the sample names to be included in the variant selection process. Ensure that these names match exactly with those in your VCF file.

gatk SelectVariants \
-V /Users/yaz/Desktop/r3/Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf \
-sn KAG-01tl -sn KAG-02tl -sn KAG-03tl -sn KAG-04tl \
-sn LIC-01tl -sn LIC-02tl -sn LIC-03tl \
-sn MOD-01tl -sn MOD-02tl -sn MOD-03tl -sn MOD-04tl \
-sn FRE-01tl -sn FRE-02tl -sn FRE-03tl -sn FRE-04tl -sn FRE-05tl -sn FRE-06tl -sn FRE-07tl -sn FRE-08tl \
-sn OCH-01tl -sn OCH-02tl -sn OCH-03tl -sn OCH-04tl -sn OCH-05tl -sn OCH-06tl -sn OCH-07tl -sn OCH-08tl \
-sn HAB-01tl -sn HAB-02tl -sn HAB-03tl \
-sn PIZ-03dl -sn PIZ-04dl -sn PIZ-05dl -sn PIZ-06dl -sn PIZ-08dl -sn PIZ-09dl -sn PIZ-11dl \
-sn VLH-01dl -sn VLH-02dl -sn VLH-03dl -sn VLH-05dl -sn VLH-06dl -sn VLH-07dl \
-sn KEH-01tl -sn KEH-04tl -sn KEH-05tl -sn KEH-06tl -sn KEH-07tl -sn KEH-08tl -sn KEH-09tl -sn KEH-10tl \
-sn BZD-01tl -sn BZD-02tl -sn BZD-03tl -sn BZD-04tl -sn BZD-05tl -sn BZD-06tl -sn BZD-07tl -sn BZD-08tl \
-O /Users/yaz/Desktop/r3/my_filtered_output.vcf 

# Deactivate environment                                 
                                                                               
conda deactivate                      
