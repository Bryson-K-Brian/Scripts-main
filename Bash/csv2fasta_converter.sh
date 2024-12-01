#!/bin/bash

# Enable tab completion for file selection
read -e -p "Enter the input CSV file name (with extension): " input_file

# Check if the input file exists
if [[ ! -f "$input_file" ]]; then
    echo "Input file '$input_file' not found!"
    exit 1
fi

# Enable tab completion for the output file name
read -e -p "Enter the output FASTA file name (with extension): " output_file

# Create or overwrite the output FASTA file
> "$output_file"

# Read the CSV file line by line
while IFS=',' read -r id sequence; do
    # Skip the header row if present
    if [[ "$id" == "id" && "$sequence" == "sequence" ]]; then
        continue
    fi
    # Add '>' to the header if it's missing
    if [[ "$id" != ">"* ]]; then
        id=">$id"
    fi
    # Write to FASTA file in the proper format
    echo "$id" >> "$output_file"
    echo "$sequence" >> "$output_file"
done < "$input_file"

echo "Conversion completed. Output saved in '$output_file'."