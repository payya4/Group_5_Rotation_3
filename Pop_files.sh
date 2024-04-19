#USING PYTHON 3.12.0

# List of your sample names, ensuring proper format for easy parsing
sample_names = [
    "KAG-01tl", "KAG-02tl",
    # Add the rest of your sample names here
    "KAG-03tl", "MOD-01tl", "FRE-01tl"
]

# Define the path for the output file
output_file_path = "individuals_population_ids.txt"

# Open a new text file to write the output
with open(output_file_path, "w") as file:
    for sample in sample_names:
        # Extract the population ID from the sample name
        population_id = sample.split("-")[0]
        # Write the sample name and population ID to the file in a tab-delimited format
        file.write(f"{sample}\t{population_id}\n")

print(f"File '{output_file_path}' written successfully with {len(sample_names)} entries.")
