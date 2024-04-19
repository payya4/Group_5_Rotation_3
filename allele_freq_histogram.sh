# Load necessary libraries
library(dplyr)
library(ggplot2)

# Read in allele frequency data from two different populations
sample_pop1 <- read.table("/path/to/freq_with.txt", header = TRUE) # Modify path as needed
sample_pop2 <- read.table("/path/to/freq_without.txt", header = TRUE) # Modify path as needed

# Merge datasets by chromosome and position
hybrid <- merge(sample_pop1, sample_pop2, by = c("CHROM", "POS"), all = TRUE) %>%
  select(c("CHROM", "POS", "AF.x", "AF.y"))

# Handle missing values and calculate mean allele frequencies
sim <- na.omit(hybrid)
sim$AF.mean <- (sim$AF.x + sim$AF.y) / 2

# Create a PDF file for the histogram output
pdf("Histogram_of_Hybrid_Simulation.pdf", height = 5, width = 9)

# Generate histogram of mean allele frequencies
ggplot(sim, aes(x = AF.mean)) +
  geom_histogram(binwidth = 0.05, fill = "blue", color = "black") +
  labs(title = "Histogram of Hybrid Simulation",
       x = "Allele Frequency",
       y = "Count") + 
  theme_minimal()

# Close the PDF device
dev.off()
