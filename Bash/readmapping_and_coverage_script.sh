#!/bin/bash

echo "### Read Mapping and Coverage Analysis ###"

# Prompt for input files with tab completion
read -e -p "Enter the path to the contigs file (FASTA): " CONTIGS
if [ ! -f "$CONTIGS" ]; then
    echo "Error: Contigs file not found!"
    exit 1
fi

read -e -p "Enter the path to the first reads file (FASTQ): " READS_R1
if [ ! -f "$READS_R1" ]; then
    echo "Error: Reads file not found!"
    exit 1
fi

read -p "Are you using paired-end reads? (yes/no): " PAIR_END
if [[ "$PAIR_END" =~ ^[Yy] ]]; then
    read -e -p "Enter the path to the second reads file (FASTQ): " READS_R2
    if [ ! -f "$READS_R2" ]; then
        echo "Error: Second reads file not found!"
        exit 1
    fi
else
    READS_R2=""
fi

read -p "Enter the sample name: " SAMPLE_NAME
read -p "Enter the number of threads to use (default: 4): " THREADS
THREADS=${THREADS:-4}

read -p "Enter the minimum coverage threshold (default: 5): " COVERAGE_THRESHOLD
COVERAGE_THRESHOLD=${COVERAGE_THRESHOLD:-5}

# Create output directory
OUTPUT_DIR="coverage/$SAMPLE_NAME"
mkdir -p "$OUTPUT_DIR"

PREFIX="$OUTPUT_DIR/mapping_output"
echo "Using $THREADS threads and filtering with coverage threshold of $COVERAGE_THRESHOLD."
echo "Output will be stored in: $OUTPUT_DIR"

# Index the contigs if not already indexed
if [ ! -f "${CONTIGS}.bwt" ]; then
    echo "Indexing contigs..."
    bwa index "$CONTIGS"
fi

# Align reads to contigs
if [ -z "$READS_R2" ]; then
    echo "Performing single-end mapping..."
    bwa mem -t "$THREADS" "$CONTIGS" "$READS_R1" | samtools view -@ "$THREADS" -bS - | samtools sort -@ "$THREADS" -o "$PREFIX.sorted.bam"
else
    echo "Performing paired-end mapping..."
    bwa mem -t "$THREADS" "$CONTIGS" "$READS_R1" "$READS_R2" | samtools view -@ "$THREADS" -bS - | samtools sort -@ "$THREADS" -o "$PREFIX.sorted.bam"
fi

# Index the sorted BAM file
samtools index "$PREFIX.sorted.bam"

# Compute depth of coverage and save as TSV
echo -e "Contig\tPosition\tCoverage" > "$PREFIX.coverage_per_position.tsv"
samtools depth -a "$PREFIX.sorted.bam" >> "$PREFIX.coverage_per_position.tsv"

# Filter for low coverage
awk -v threshold=$COVERAGE_THRESHOLD '$3 >= threshold' "$PREFIX.coverage_per_position.tsv" > "$PREFIX.filtered_coverage.tsv"

# Compute average coverage per contig and save to TSV
echo -e "Contig\tAverage_Coverage" > "$PREFIX.avg_coverage_per_contig.tsv"
awk '{sum[$1]+=$3; count[$1]++} END {for (c in sum) print c "\t" sum[c]/count[c]}' "$PREFIX.coverage_per_position.tsv" >> "$PREFIX.avg_coverage_per_contig.tsv"

# Compute overall average coverage
OVERALL_AVG=$(awk '{sum+=$3; count++} END {if (count > 0) print sum/count; else print "N/A"}' "$PREFIX.coverage_per_position.tsv")
echo -e "Overall_Average_Coverage\t$OVERALL_AVG" > "$PREFIX.overall_avg_coverage.tsv"

# Generate coverage visualization using Python
python3 - <<EOF
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Load depth data
depth_data = pd.read_csv("$PREFIX.coverage_per_position.tsv", sep="\t", header=0)

# Plot coverage depth across contigs
plt.figure(figsize=(10, 6))
for contig in depth_data["Contig"].unique():
    contig_data = depth_data[depth_data["Contig"] == contig]
    plt.plot(contig_data["Position"], contig_data["Coverage"], label=f"Contig {contig}")

plt.xlabel("Position")
plt.ylabel("Coverage Depth")
plt.title("Coverage Depth Across Contigs")
plt.legend()
plt.tight_layout()
plt.savefig("$PREFIX.coverage_plot.png")
plt.close()

# Generate histogram of coverage distribution
plt.figure(figsize=(8, 5))
sns.histplot(depth_data["Coverage"], bins=50, kde=True, color="skyblue")
plt.xlabel("Coverage Depth")
plt.ylabel("Frequency")
plt.title("Histogram of Coverage Depth Distribution")
plt.savefig("$PREFIX.coverage_histogram.png")
plt.close()

EOF

echo "Mapping, filtering, and depth analysis complete."
echo "Coverage per position saved in: $PREFIX.coverage_per_position.tsv"
echo "Filtered coverage (>$COVERAGE_THRESHOLD) saved in: $PREFIX.filtered_coverage.tsv"
echo "Average coverage per contig saved in: $PREFIX.avg_coverage_per_contig.tsv"
echo "Overall average coverage saved in: $PREFIX.overall_avg_coverage.tsv"
echo "Coverage depth plot saved as: $PREFIX.coverage_plot.png"
echo "Coverage histogram saved as: $PREFIX.coverage_histogram.png"
echo "Results are stored in: $OUTPUT_DIR"