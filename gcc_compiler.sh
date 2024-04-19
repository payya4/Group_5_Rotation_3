## Usage of poly_freq.c for Allele Frequency Calculations
###The `poly_freq.c` script is used to calculate allele frequencies from a VCF file. Below are the steps to compile and run the script, along with examples of how to use it both with and without considering outlier data.

##Compiling the Script
##First, compile the script using GCC. Navigate to the directory containing `poly_freq.c` and run:
gcc /Users/yaz/Desktop/r3/poly_freq.c -o poly_freq -lm


#Calculating Allele Frequencies With Outliers:
./poly_freq -vcf /Users/yaz/Desktop/r3/my_filtered_output.vcf -pops /Users/yaz/Desktop/r3/with_individuals_population_ids.txt  > /Users/yaz/Desktop/r3/freq_with


#Calculating Allele Frequencies Without Outliers:
./poly_freq -vcf /Users/yaz/Desktop/r3/my_filtered_output_no_outliers.vcf -pops /Users/yaz/Desktop/r3/without_individuals_population_ids_updated.txt > /Users/yaz/Desktop/r3/freq_without


