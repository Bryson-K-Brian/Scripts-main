#!/bin/bash

# Ensure required tools are available
for cmd in seqtk parallel; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: $cmd is required but not installed. Install it with 'conda install -c bioconda $cmd' or 'brew install $cmd'."
        exit 1
    fi
done

# Enable tab completion for user input
shopt -s progcomp

# Prompt for input directories with tab completion
read -e -p "Enter the directory containing Kraken2 output files: " KRAKEN_OUT_DIR
read -e -p "Enter the directory containing original contigs: " CONTIGS_DIR
read -e -p "Enter the output directory for decontaminated contigs: " CLEAN_DIR

# Create output directory if it does not exist
mkdir -p "$CLEAN_DIR"

# Prompt for species name
read -p "Enter the species name (e.g., 'Aeromonas veronii'): " SPECIES

# Get available CPU cores and set parallel jobs
CPU_CORES=$(nproc)
JOBS=$((CPU_CORES / 2))  # Adjust parallel jobs (use half the cores)

export SPECIES KRAKEN_OUT_DIR CONTIGS_DIR CLEAN_DIR

# Function to process a single Kraken2 report
process_report() {
    REPORT="$1"
    BASE=$(basename "$REPORT" "_output.tsv")
    CONTIG_FILE="$CONTIGS_DIR/$BASE.fasta"
    
    if [[ ! -f "$CONTIG_FILE" ]]; then
        echo "Warning: Contigs file $CONTIG_FILE not found, skipping..."
        return
    fi

    # Extract contig IDs classified at the species level (including subspecies and strains)
    grep -i -E "^$SPECIES(\s|$)" "$REPORT" | cut -f2 > "${BASE}_target_contigs.list"

    # Filter the contigs
    if [[ -s "${BASE}_target_contigs.list" ]]; then
        seqtk subseq "$CONTIG_FILE" "${BASE}_target_contigs.list" > "$CLEAN_DIR/${BASE}_decontaminated.fasta"
        echo "Filtered contigs for $BASE saved to $CLEAN_DIR/${BASE}_decontaminated.fasta"
    else
        echo "No contigs matched for $BASE, skipping..."
    fi

    # Clean up
    rm -f "${BASE}_target_contigs.list"
}

export -f process_report

# Parallel processing of Kraken2 reports
find "$KRAKEN_OUT_DIR" -name "*_output.tsv" | parallel -j "$JOBS" process_report {}

echo "Parallel decontamination process completed."