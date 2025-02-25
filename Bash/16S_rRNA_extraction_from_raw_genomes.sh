#!/bin/bash

# Enable tab completion for input
read -e -p "Enter the directory containing genome FASTA files: " INPUT_DIR
read -e -p "Enter the output file name: " OUTPUT_FILE

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Directory '$INPUT_DIR' does not exist."
    exit 1
fi

# Create or clear the output file
> "$OUTPUT_FILE"

# Loop through all FASTA files in the directory
for genome in "$INPUT_DIR"/*.fasta; do
    # Extract genome base name (without path and extension)
    genome_name=$(basename "$genome" .fasta)

    echo "Processing $genome_name..."

    # Run Barrnap to predict rRNA genes
    barrnap --kingdom bacteria "$genome" > "${genome_name}_rRNA.gff"

    # Check if the GFF file contains 16S rRNA annotations
    if grep -q "16S" "${genome_name}_rRNA.gff"; then
        # Extract only the first 16S rRNA gene entry
        grep "16S" "${genome_name}_rRNA.gff" | awk '!seen[$1]++' > "${genome_name}_16S.gff"

        # Extract the 16S sequence using bedtools
        bedtools getfasta -fi "$genome" -bed "${genome_name}_16S.gff" -fo "${genome_name}_16S.fasta"

        # Rename header to just the genome name
        awk -v name="$genome_name" '/^>/ {print ">" name; next} {print}' "${genome_name}_16S.fasta" >> "$OUTPUT_FILE"
    else
        echo "Warning: No 16S rRNA gene found in $genome_name"
    fi

    # Cleanup temporary files
    rm -f "${genome_name}_rRNA.gff" "${genome_name}_16S.gff" "${genome_name}_16S.fasta"
done

echo "✅ 16S rRNA sequences saved to '$OUTPUT_FILE'"