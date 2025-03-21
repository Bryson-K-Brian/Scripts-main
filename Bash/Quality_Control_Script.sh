#!/bin/bash

# Prompt for the input directory containing raw reads
read -e -p "Enter the input directory containing raw reads: " input_dir

# Check if the input directory exists
if [[ ! -d "$input_dir" ]]; then
    echo "Input directory '$input_dir' does not exist. Exiting."
    exit 1
fi

# Prompt for the output directory prefix
read -e -p "Enter the output directory (e.g., 'trimmed_results'): " output_directory

# Create the output directory
output_dir="$output_directory"
mkdir -p "$output_dir"

# Check if there are matching FASTQ files
if ! ls "$input_dir"/*01.fastq.gz &>/dev/null; then
    echo "No matching FASTQ files found in '$input_dir'. Exiting."
    exit 1
fi

# Loop through Read1 files and process them
for f in "$input_dir"/*R1_001.fastq.gz; do
    # Define Read1 and Read2
    RD1="$f"
    RD2="${f%R1_001.fastq.gz}R2_001.fastq.gz"

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

    # Define repair output files
    rpr1="${output_dir}/$(basename "${f%R1_001.fastq.gz}repaired_R1_001.fastq.gz")"
    rpr2="${output_dir}/$(basename "${f%R1_001.fastq.gz}repaired_R2_001.fastq.gz")"

    # Repair read pairs for disorder and missing mates
    echo "Running Repair for sample $(basename "$f" | sed 's/_S.*//')..."
    repair.sh in="$RD1" in2="$RD2" out="$rpr1" out2="$rpr2"
    
    # Define trimmed output files
    trm1="${output_dir}/$(basename "${f%R1_001.fastq.gz}trimmed_R1_001.fastq.gz")"
    trm2="${output_dir}/$(basename "${f%R1_001.fastq.gz}trimmed_R2_001.fastq.gz")"

    # Quality Control using BBDUK
    echo "Running BBDUK for sample $(basename "$f" | sed 's/_S.*//')..."
    bbduk.sh in="$rpr1" in2="$rpr2" out="$trm1" out2="$trm2" trimq=20 qtrim=rl minlength=50 minbasequality=0 ref="adapters"

    # Log completion for the sample
    echo "Sample $(basename "$f" | sed 's/_S.*//') Complete."
    echo ""
done

echo "All trimmed files have been saved to the '$output_dir' directory."