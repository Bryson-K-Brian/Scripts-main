from Bio import SeqIO
from collections import defaultdict

input_file = "BVBRC_genome_sequence.fasta"
output_dir = "output/"

# Dictionary to hold sequences by genome
genomes = defaultdict(list)

# Parse the FASTA file
for record in SeqIO.parse(input_file, "fasta"):
    # Extract the content within square brackets, split by '|' and get the first part
    if "[" in record.description:
        genome_id = record.description.split("[")[1].split("|")[0].strip()
    else:
        genome_id = "unknown"

    genomes[genome_id].append(record)

# Write each genome's contigs to a separate file
for genome_id, records in genomes.items():
    sanitized_genome_id = genome_id.replace(' ', '_')  # Replace spaces with underscores
    output_file = f"{output_dir}{sanitized_genome_id}.fasta"
    with open(output_file, "w") as f:
        SeqIO.write(records, f, "fasta")
