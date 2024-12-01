#!/bin/bash

# Prompt for the input directory containing reads/contigs
read -e -p "Enter the input directory containing the reads or contigs: " input_dir

# Check if the input directory exists
if [[ ! -d "$input_dir" ]]; then
    echo "Input directory '$input_dir' does not exist. Exiting."
    exit 1
fi

# Prompt for the Kraken2 database selection
echo "Select the Kraken2 database:"
select db in "viral_db" "standard_db" "other_db"; do
    case $db in
        viral_db)
            echo "You selected the viral database"
            break
            ;;
        standard_db)
            echo "You selected the standard database"
            break
            ;;
        other_db)
            echo "You selected another database (make sure it is correctly installed)"
            break
            ;;
        *)
            echo "Invalid selection. Please choose a valid database."
            ;;
    esac
done

# Prompt for the type of input (paired-end, single-end, or contigs)
echo "Select the type of input data:"
select read_type in "paired_end" "single_end" "contigs"; do
    case $read_type in
        paired_end)
            echo "You selected paired-end reads"
            break
            ;;
        single_end)
            echo "You selected single-end reads"
            break
            ;;
        contigs)
            echo "You selected contigs"
            break
            ;;
        *)
            echo "Invalid selection. Please choose a valid read type."
            ;;
    esac
done

# Create the output directory for Kraken2 results
mkdir -p taxid

# Navigate to the input directory
cd "$input_dir" || exit 1

# Handle paired-end reads
if [[ "$read_type" == "paired_end" ]]; then
    # Check for paired FASTQ files
    if ! ls *R1_001.fastq.gz &>/dev/null; then
        echo "No matching FASTQ files (*R1_001.fastq.gz) found in '$input_dir'. Exiting."
        exit 1
    fi

    # Process each Read1 file
    for f in *R1_001.fastq.gz; do
        # Define Read1 and Read2
        RD1="$f"
        RD2="${f%_R1_001.fastq.gz}_R2_001.fastq.gz"

        # Log the file pairs
        echo ""
        echo "Processing paired-end sample: $(basename "$f" | sed 's/_R1_001.fastq.gz//')"
        echo "Read1: $RD1"
        echo "Read2: $RD2"
        echo ""

        # Check if both Read1 and Read2 exist
        if [[ ! -f "$RD1" || ! -f "$RD2" ]]; then
            echo "Error: Missing pair for $RD1 and $RD2. Skipping."
            continue
        fi

        # Run Kraken2 for paired-end taxonomic identification
        echo "Running Kraken2 for paired-end sample $(basename "$f" | sed 's/_R1_001.fastq.gz//') using the $db database..."
        kraken2 --db "$db" --threads 20 --output "${f%_S*}.tsv" --report --minimum-base-quality 20 --paired --use-names "$RD1" "$RD2"

        # Move the results to the 'taxid' folder
        mv *.tsv taxid

        echo ""
        echo "Sample $(basename "$f" | sed 's/_R1_001.fastq.gz//') Kraken2 taxonomic identification complete."
        echo ""
    done

# Handle single-end reads
elif [[ "$read_type" == "single_end" ]]; then
    # Check for single-end FASTQ files
    if ! ls *.fastq.gz &>/dev/null; then
        echo "No matching FASTQ files (*.fastq.gz) found in '$input_dir'. Exiting."
        exit 1
    fi

    # Process each single-end FASTQ file
    for f in *.fastq.gz; do
        # Log the file
        echo ""
        echo "Processing single-end sample: $(basename "$f" | sed 's/.fastq.gz//')"
        echo "Read: $f"
        echo ""

        # Run Kraken2 for single-end taxonomic identification
        echo "Running Kraken2 for single-end sample $(basename "$f" | sed 's/.fastq.gz//') using the $db database..."
        kraken2 --db "$db" --threads 20 --output "${f%_S*}.tsv" --report --minimum-base-quality 20 --use-names "$f"

        # Move the results to the 'taxid' folder
        mv *.tsv taxid

        echo ""
        echo "Sample $(basename "$f" | sed 's/.fastq.gz//') Kraken2 taxonomic identification complete."
        echo ""
    done

# Handle contig files
elif [[ "$read_type" == "contigs" ]]; then
    # Check for contig files (assumes they are in FASTA format)
    if ! ls *.fasta &>/dev/null; then
        echo "No matching contig files (*.fasta) found in '$input_dir'. Exiting."
        exit 1
    fi

    # Process each contig file
    for f in *.fasta; do
        # Log the file
        echo ""
        echo "Processing contig sample: $(basename "$f" | sed 's/.fasta//')"
        echo "Contig file: $f"
        echo ""

        # Run Kraken2 for contig taxonomic identification
        echo "Running Kraken2 for contig sample $(basename "$f" | sed 's/.fasta//') using the $db database..."
        kraken2 --db "$db" --threads 20 --output "${f%_S*}.tsv" --report --minimum-base-quality 20 --use-names "$f"

        # Move the results to the 'taxid' folder
        mv *.tsv taxid

        echo ""
        echo "Sample $(basename "$f" | sed 's/.fasta//') Kraken2 taxonomic identification complete."
        echo ""
    done
fi

# Navigate back to the original directory
cd - >/dev/null || exit 1

echo "All samples processed. Taxonomic results are saved in the 'taxid' directory."