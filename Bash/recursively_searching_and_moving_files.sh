#!/bin/bash

# Prompt for the directory to search
read -e -p "Enter the directory to search in: " search_dir

# Check if the directory exists
if [[ ! -d "$search_dir" ]]; then
    echo "Directory '$search_dir' does not exist. Exiting."
    exit 1
fi

# Prompt for the file pattern (e.g., *.fastq.gz, *.bam, *.txt)
read -p "Enter the file pattern to search for (e.g., *.fastq.gz): " file_pattern

# Prompt for the destination directory
read -e -p "Enter the destination directory: " dest_dir

# Create the destination directory if it doesn't exist
mkdir -p "$dest_dir"

# Find and move matching files
find "$search_dir" -type f -name "$file_pattern" -exec mv {} "$dest_dir" \;

echo "All matching files have been moved to '$dest_dir'."