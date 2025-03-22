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

read -p "Enter the amount of RAM to allocate (in GB): " SYSTEM_RAM
SYSTEM_RAM=${SYSTEM_RAM:-8}  # Default to 8GB if no input provided
ALLOCATED_RAM=$((SYSTEM_RAM * 1024))  # Convert to MB

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

# Map the user selection to Unicycler mode
case $assembly_type in
    isolate|plasmid)
        unicycler_mode="normal"
        ;;
    meta|metaviral|metaplasmid)
        unicycler_mode="bold"
        ;;
    rna|rna_viral)
        unicycler_mode="conservative"
        ;;
    *)
        echo "Invalid assembly type for Unicycler. Defaulting to 'normal'."
        unicycler_mode="normal"
        ;;
esac

# Prompt for the output directory
read -e -p "Enter the output directory prefix (e.g., wgs_assembly): " output_dir

# Create the output directory
mkdir -p "$output_dir"

# Navigate to the input directory
cd "$input_dir" || exit 1

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
    echo "Processing sample: $(basename "$f" | sed 's/_R1_001.fastq.gz//')"
    echo "Read1: $RD1"
    echo "Read2: $RD2"
    echo ""

    # Check if both Read1 and Read2 exist
    if [[ ! -f "$RD1" || ! -f "$RD2" ]]; then
        echo "Error: Missing pair for $RD1 and $RD2. Skipping."
        continue
    fi

    # Define output folder for this sample
    sample_name="${f%_R1_001.fastq.gz}"
    sample_output_dir="$output_dir/${sample_name}"

    # Create output directory for the sample
    mkdir -p "$sample_output_dir"

    # Run assembler (SPAdes or Unicycler)
    if [[ "$assembler" == "spades" ]]; then
        echo "Running SPAdes for sample $sample_name with assembly type $assembly_type"
        spades.py -1 "$RD1" -2 "$RD2" -o "${sample_name}_assembly/spades_output" --threads "$THREADS" --$assembly_type 
        # Rename SPAdes contigs file
        mv ${sample_name}_assembly/spades_output/contigs.fasta ${sample_name}_assembly/${sample_name}_spades.fasta
    elif [[ "$assembler" == "unicycler" ]]; then
        echo "Running Unicycler for sample $sample_name"
        unicycler --verbosity 1 -1 "$RD1" -2 "$RD2" -o "$sample_output_dir/unicycler_output" --threads "$THREADS" --mode "$unicycler_mode" --spades_options "-m $ALLOCATED_RAM"
        if [[ -f "$sample_output_dir/unicycler_output/assembly.fasta" ]]; then
            mv "$sample_output_dir/unicycler_output/assembly.fasta" "$sample_output_dir/${sample_name}_unicycler.fasta"
        else
            echo "Error: Unicycler did not generate an assembly.fasta file for $sample_name"
        fi
    fi

    echo ""
    echo "Sample $sample_name processing complete."
    echo ""
done

# Navigate back to the original directory
cd - >/dev/null || exit 1

echo "All samples processed. Results are saved in the '$output_dir' directory."