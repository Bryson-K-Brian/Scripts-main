#!/bin/bash

echo "### Genome Coverage & Sequence Identity Analysis ###"

# Prompt for input files with tab completion
read -e -p "Enter the reference genome (FASTA): " REF_GENOME
if [ ! -f "$REF_GENOME" ]; then
    echo "Error: Reference genome not found!"
    exit 1
fi

read -e -p "Enter the reads file (FASTQ): " READS
if [ ! -f "$READS" ]; then
    echo "Error: Reads file not found!"
    exit 1
fi

read -p "Are the reads paired-end? (yes/no): " PAIR_END
if [[ "$PAIR_END" =~ ^[Yy] ]]; then
    read -e -p "Enter the second reads file (FASTQ): " READS2
    if [ ! -f "$READS2" ]; then
        echo "Error: Second reads file not found!"
        exit 1
    fi
else
    READS2=""
fi

read -p "Enter the sample name: " SAMPLE_NAME
read -p "Enter the number of threads (default: 4): " THREADS
THREADS=${THREADS:-4}

# Create output directory
OUTPUT_DIR="coverage_identity/$SAMPLE_NAME"
mkdir -p "$OUTPUT_DIR"

PREFIX="$OUTPUT_DIR/${SAMPLE_NAME}"
echo "Using $THREADS threads..."
echo "Results will be stored in: $OUTPUT_DIR"

# Index the reference genome if not indexed
if [ ! -f "${REF_GENOME}.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index "$REF_GENOME"
fi

# Align reads to reference
if [ -z "$READS2" ]; then
    echo "Performing single-end mapping..."
    bwa mem -t "$THREADS" "$REF_GENOME" "$READS" | samtools view -@ "$THREADS" -bS - | samtools sort -@ "$THREADS" -o "$PREFIX.sorted.bam"
else
    echo "Performing paired-end mapping..."
    bwa mem -t "$THREADS" "$REF_GENOME" "$READS" "$READS2" | samtools view -@ "$THREADS" -bS - | samtools sort -@ "$THREADS" -o "$PREFIX.sorted.bam"
fi

# Index the sorted BAM file
samtools index "$PREFIX.sorted.bam"

# Compute coverage per base
samtools depth -a "$PREFIX.sorted.bam" > "$PREFIX.depth.txt"

# Calculate genome coverage
TOTAL_COVERED=$(awk '$3>0 {covered++} END {print covered+0}' "$PREFIX.depth.txt")
GENOME_SIZE=$(awk '{sum+=$2} END {print (sum==0)?1:sum}' <(samtools faidx "$REF_GENOME"))
TOTAL_GENOME_BASES=$(awk '{sum+=$3} END {print (sum==0)?1:sum}' "$PREFIX.depth.txt")

COVERAGE_PERCENT=$(echo "scale=2; ($TOTAL_COVERED/$GENOME_SIZE)*100" | bc || echo "0")
REF_GENOME_COVERAGE=$(echo "scale=2; ($TOTAL_GENOME_BASES/$GENOME_SIZE)*100" | bc || echo "0")
echo "Reference Genome Coverage: $REF_GENOME_COVERAGE%" | tee -a "$PREFIX.coverage_summary.txt"

# Run Mummer dnadiff for sequence identity
echo "Calculating sequence identity..."
dnadiff -p "$PREFIX.dnadiff" "$REF_GENOME" "$READS"

# Extract identity from dnadiff report
IDENTITY=$(grep "AvgIdentity" "$PREFIX.dnadiff.report" | awk '{print $3}')
echo "Sequence Identity: $IDENTITY%" | tee -a "$PREFIX.coverage_summary.txt"

echo "Analysis complete!"
echo "Results saved in: $OUTPUT_DIR"