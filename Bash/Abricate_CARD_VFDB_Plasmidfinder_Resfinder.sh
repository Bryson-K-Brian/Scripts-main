#!/bin/bash

# Prompt for the input folder with tab completion
read -e -p "Enter the input folder containing FASTA files: " input_folder

# Check if the input folder exists
if [[ ! -d "$input_folder" ]]; then
    echo "Input folder '$input_folder' does not exist. Exiting."
    exit 1
fi

# Prompt for the output prefix
read -e -p "Enter the output folder prefix (e.g., 'abricate_results'): " output_prefix

# Create the output directory
output_dir="${output_prefix}_abricate"
mkdir -p "$output_dir"

# List of databases to process
databases=("card" "plasmidfinder" "vfdb" "resfinder")

# Check if there are any FASTA files in the input folder
if ! ls "$input_folder"/*.fasta &>/dev/null; then
    echo "No FASTA files found in the folder '$input_folder'. Exiting."
    exit 1
fi

# Function to run abricate for a given file and database
run_abricate() {
    local file="$1"
    local db="$2"
    local base_name
    base_name=$(basename "${file%.fasta}")
    local output_file="$output_dir/${base_name}_${db}.csv"

    echo "Processing $file with database $db..."
    abricate --csv --db "$db" "$file" > "$output_file"
    echo "Sample $base_name with database $db done."
}

export -f run_abricate # Export the function for parallel use

# Loop through databases and process files
for db in "${databases[@]}"; do
    echo "Starting abricate analysis for database: $db"

    # Use GNU Parallel if available, fallback to sequential processing
    if command -v parallel &>/dev/null; then
        find "$input_folder" -name "*.fasta" | parallel run_abricate {} "$db"
    else
        for f in "$input_folder"/*.fasta; do
            run_abricate "$f" "$db"
        done
    fi

    echo "Completed abricate analysis for database: $db"
done

echo "All analyses completed. Results are in the '$output_dir' directory."