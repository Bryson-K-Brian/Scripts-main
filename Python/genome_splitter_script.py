from Bio import SeqIO
from collections import defaultdict
import os
import readline
import glob
import re

# Function for tab completion
def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]

readline.set_completer(complete)
readline.parse_and_bind("tab: complete")

# Ask for input file
input_file = input("Enter the path to the input FASTA file: ")
while not os.path.isfile(input_file):
    print("File not found. Please enter a valid path.")
    input_file = input("Enter the path to the input FASTA file: ")

output_dir = "output/"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Dictionary to hold sequences by genome
genomes = defaultdict(list)

# Parse the FASTA file
for record in SeqIO.parse(input_file, "fasta"):
    # Extract genome identifier from description
    if "[" in record.description:
        genome_id = record.description.split("[")[1].split("|")[0].strip()
    else:
        genome_id = "unknown"
    
    # Sanitize genome ID to be filesystem safe
    sanitized_genome_id = re.sub(r'[^a-zA-Z0-9._-]', '_', genome_id)
    
    genomes[sanitized_genome_id].append(record)

# Write each genome's sequences to a separate FASTA file
for genome_id, records in genomes.items():
    output_file = os.path.join(output_dir, f"{genome_id}.fasta")
    
    # Ensure the directory exists before writing the file
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    with open(output_file, "w") as f:
        SeqIO.write(records, f, "fasta")

print(f"Processing complete. Separated genomes are saved in '{output_dir}'")
