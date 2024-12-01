#!/bin/bash

# Prompt for the input directory containing FASTQ files
read -e -p "Enter the input directory containing FASTQ files: " input_dir

# Check if the input directory exists
if [[ ! -d "$input_dir" ]]; then
    echo "Input directory '$input_dir' does not exist. Exiting."
    exit 1
fi

# Prompt for the reference genome file
read -e -p "Enter the reference genome file (e.g., *.fna): " reference_genome

# Check if the reference genome file exists
if [[ ! -f "$reference_genome" ]]; then
    echo "Reference genome file '$reference_genome' does not exist. Exiting."
    exit 1
fi

# Create the output directory
output_dir="readmapping"
mkdir -p "$output_dir"

# Navigate to the input directory
cd "$input_dir" || exit 1

# Check for paired FASTQ files
if ! ls *_R1_001.fastq.gz &>/dev/null; then
    echo "No matching FASTQ files (*_R1_001.fastq.gz) found in '$input_dir'. Exiting."
    exit 1
fi

# Process each Read1 file
for f in *_R1_001.fastq.gz; do
    # Define Read1 and Read2
    RD1="$f"
    RD2="${f%_R1_001.fastq.gz}_R2_001.fastq.gz"

    # Log the file pairs
    echo ""
    echo "Processing sample: $(basename "$f" | sed 's/_R1_001.fastq.gz//')"
    echo "Read1: $RD1"
    echo "Read2: $RD2"
    echo ""

    # Check if both Read1 and Read2 exist
    if [[ ! -f "$RD1" || ! -f "$RD2" ]]; then
        echo "Error: Missing pair for $RD1 and $RD2. Skipping."
        continue
    fi

    # Define output filenames
    sample_name="${f%_R1_001.fastq.gz}"
    sam_file="${sample_name}.sam"
    bam_file="${sample_name}.bam"
    bam_12F="${sample_name}_12F.bam"
    sorted_bam="${sample_name}_12F_sorted.bam"
    depth_file="${sample_name}_12F_sorted.bam.txt"

    # Perform read mapping and processing
    echo "Running BWA-MEM for sample $sample_name..."
    bwa mem -t 6 "$reference_genome" "$RD1" "$RD2" -o "$sam_file"

    echo "Converting SAM to BAM..."
    samtools view -b "$sam_file" > "$bam_file"

    echo "Filtering BAM (flag 12)..."
    samtools view -bF 12 "$bam_file" > "$bam_12F"

    echo "Sorting BAM..."
    samtools sort "$bam_12F" > "$sorted_bam"

    echo "Generating sequencing depth file..."
    samtools depth "$sorted_bam" > "$depth_file"

    # Organize output
    echo "Organizing output files..."
    sample_dir="$output_dir/$sample_name"
    mkdir -p "$sample_dir"
    mv "$sample_name"* "$sample_dir"

    echo "Sample $sample_name processing complete."
    echo ""
done

# Navigate back to the original directory
cd - >/dev/null || exit 1

echo "All samples processed. Results are saved in the '$output_dir' directory."