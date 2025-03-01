#!/bin/bash

# Prompt for the input directory containing trimmed reads
read -e -p "Enter the input directory containing trimmed reads: " input_dir

# Check if the input directory exists
if [[ ! -d "$input_dir" ]]; then
    echo "Input directory '$input_dir' does not exist. Exiting."
    exit 1
fi

#Prompt for processing inputs
read -p "Enter the number of threads to use: " THREADS

# Prompt for the assembler selection (SPAdes or Unicycler)
echo "Select the assembler:"
select assembler in "spades" "unicycler"; do
    case $assembler in
        spades)
            echo "You selected SPAdes"
            break
            ;;
        unicycler)
            echo "You selected Unicycler"
            break
            ;;
        *)
            echo "Invalid selection. Please choose 1 for SPAdes or 2 for Unicycler."
            ;;
    esac
done

# Prompt for the assembly type selection (plasmid, rna viral, etc.)
echo "Select the assembly type:"
select assembly_type in "rna_viral" "plasmid" "isolate" "meta" "rna" "metaviral" "metaplasmid"; do
    case $assembly_type in
        rna_viral)
            echo "You selected RNA viral assembly"
            break
            ;;
        plasmid)
            echo "You selected Plasmid assembly"
            break
            ;;
        isolate)
            echo "You selected isolate assembly"
            break
            ;;
        meta)
            echo "You selected metagenomic assembly"
            break
            ;;
        rna)
            echo "You selected rna-seq assembly"
            break
            ;;
        metaviral)
            echo "You selected isolate assembly"
            break
            ;;
        metaplasmid)
            echo "You selected isolate assembly"
            break
            ;;
        *)
            echo "Invalid selection. Please choose a valid assembly type."
            ;;
    esac
done

# Prompt for the output directory
read -e -p "Enter the output directory prefix (e.g., wgs_assembly): " output_dir

# Create the output directory
mkdir -p "$output_dir"

# Navigate to the input directory
cd "$input_dir" || exit 1

# Check for paired FASTQ files
if ! ls *R1_001_trimmed.fastq.gz &>/dev/null; then
    echo "No matching FASTQ files (*R1_001_trimmed.fastq.gz) found in '$input_dir'. Exiting."
    exit 1
fi

# Process each Read1 file
for f in *R1_001_trimmed.fastq.gz; do
    # Define Read1 and Read2
    RD1="$f"
    RD2="${f%_R1_001_trimmed.fastq.gz}_R2_001_trimmed.fastq.gz"

    # Log the file pairs
    echo ""
    echo "Processing sample: $(basename "$f" | sed 's/_R1_001_trimmed.fastq.gz//')"
    echo "Read1: $RD1"
    echo "Read2: $RD2"
    echo ""

    # Check if both Read1 and Read2 exist
    if [[ ! -f "$RD1" || ! -f "$RD2" ]]; then
        echo "Error: Missing pair for $RD1 and $RD2. Skipping."
        continue
    fi

    # Define output folder for this sample
    sample_name="${f%_R1_001_trimmed.fastq.gz}"
    sample_output_dir="$output_dir/${sample_name}"

    # Create output directory for the sample
    mkdir -p "$sample_output_dir"

    # Run assembler (SPAdes or Unicycler)
    if [[ "$assembler" == "spades" ]]; then
        echo "Running SPAdes for sample $sample_name with assembly type $assembly_type"
        spades.py -1 "$R1" -2 "$R2" -o ${sample_name}_assembly/spades_output --threads "$THREADS"
        # Rename SPAdes contigs file
        mv ${sample_name}_assembly/spades_output/contigs.fasta ${sample_name}_assembly/${sample_name}_spades.fasta
    elif [[ "$assembler" == "unicycler" ]]; then
        echo "Running Unicycler for sample $sample_name"
        unicycler --verbosity 1 -1 "$RD1" -2 "$RD2" -o ${sample_name}_assembly/unicycler_output --threads "$THREADS"
        mv ${sample_name}_assembly/unicycler_output/assembly.fasta ${sample_name}_assembly/${sample_name}_unicycler.fasta
    fi

    echo ""
    echo "Sample $sample_name processing complete."
    echo ""
done

# Navigate back to the original directory
cd - >/dev/null || exit 1

echo "All samples processed. Results are saved in the '$output_dir' directory."