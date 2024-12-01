#!/bin/bash

# Prompt for the input directory containing trimmed reads
read -e -p "Enter the input directory containing trimmed reads: " input_dir

# Check if the input directory exists
if [[ ! -d "$input_dir" ]]; then
    echo "Input directory '$input_dir' does not exist. Exiting."
    exit 1
fi

# Prompt for the output directory prefix
read -e -p "Enter the output directory prefix (e.g., 'metagenomic_results'): " output_prefix

# Create the output directory
output_dir="${output_prefix}_assembly"
mkdir -p "$output_dir"

# Navigate to the input directory
cd "$input_dir" || exit 1

# Check if there are matching FASTQ files
if ! ls *R1_001_trimmed.fastq.gz &>/dev/null; then
    echo "No trimmed FASTQ files found in '$input_dir'. Exiting."
    exit 1
fi

# Loop through Read1 files and process them
for f in *R1_001_trimmed.fastq.gz; do
    # Define Read1 and Read2
    RD1="$f"
    RD2="${f%R1_001_trimmed.fastq.gz}R2_001_trimmed.fastq.gz"

    # Log Read1 and Read2
    echo ""
    echo "Read1 is $RD1"
    echo "Read2 is $RD2"
    echo ""

    # Check if both Read1 and Read2 exist
    if [[ ! -f "$RD1" || ! -f "$RD2" ]]; then
        echo "Error: Missing pair for $RD1 and $RD2. Skipping."
        continue
    fi

    # Run Metagenomic Assembly with Megahit
    sample_name="${f%_S*}"
    output_path="../$output_dir/$sample_name"
    echo "Starting assembly for sample $sample_name..."
    megahit -1 "$RD1" -2 "$RD2" -o "$output_path"

    # Log completion for the sample
    echo "Sample $sample_name Complete."
    echo ""
done

# Navigate back to the original directory
cd - >/dev/null || exit 1

echo "All samples processed. Results are in the '$output_dir' directory."