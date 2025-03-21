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
output_dir="decontaminated_reads"
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
    bam_file="${sample_name}.bam"
    mapped_bam="${sample_name}_mapped.bam"
    sorted_bam="${sample_name}_mapped_sorted.bam"
    fastq_R1_out="${sample_name}_decon_R1_001.fastq.gz"
    fastq_R2_out="${sample_name}_decon_R2_001.fastq.gz"

    # Perform read mapping and processing
    echo "Running BBMap for sample $sample_name..."
    bbmap.sh ref="$reference_genome" in="$RD1" in2="$RD2" out="$bam_file" bamscript=bs.sh; bash bs.sh

    echo "Filtering BAM to retain only mapped reads..."
    samtools view -bF 4 "$bam_file" > "$mapped_bam"

    echo "Sorting BAM..."
    samtools sort "$mapped_bam" > "$sorted_bam"

    echo "Extracting mapped reads to FASTQ..."
    samtools fastq -1 "$fastq_R1_out" -2 "$fastq_R2_out" "$sorted_bam"

    # Clean up intermediate files
    echo "Cleaning up intermediate files..."
    rm "$bam_file" "$mapped_bam" "$sorted_bam"

    # Organize output
    echo "Organizing output files..."
    sample_dir="../$output_dir/$sample_name"
    mkdir -p "$sample_dir"
    mv $fastq_R1_out $fastq_R2_out "$sample_dir"

    echo "Sample $sample_name processing complete."
    echo ""
done

# Navigate back to the original directory
cd - >/dev/null || exit 1

echo "All samples processed. Decontaminated reads are saved in the '$output_dir' directory."