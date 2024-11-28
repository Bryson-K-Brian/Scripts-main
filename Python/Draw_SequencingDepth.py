import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

# Read data
input_data = sys.argv[1]
output_data = sys.argv[2]
column_names = ['Genome name', 'Genomic location', 'Depth of coverage']
data = pd.read_csv(input_data, sep='\t', header=None, names=column_names)

# retrieve data
x = data['Genomic location']
y = data['Depth of coverage']

# Calculate average depth of coverage, maximum depth, and minimum depth
average_coverage_depth = y.mean()
maximal_depth = y.max()
minimal_depth = y.min()
total_genome_length = x.max()

# Create a Graph
plt.figure(figsize=(12, 6))
plt.plot(x, y, color='green', linewidth=1, linestyle='-')

# Changing Palette
plt.xlabel('Genomic Position (Mbps)', fontsize=14)
plt.ylabel('Sequencing Depth (×)', fontsize=14)

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

plt.grid(True, linestyle='--', alpha=0.5)

# Add text
total_genome_length_str = f'{total_genome_length:,} bp'
text_str_1 = f'(1) Total genome length = {total_genome_length_str}'
text_str_3 = f'(3) Maximal depth = {maximal_depth} ×'
text_str_2 = f'(2) Average depth = {average_coverage_depth:.2f} ×'
text_str_4 = f'(4) Minimal depth = {minimal_depth} ×'

plt.gca().annotate(text_str_1, fontsize=12, xy=(0, -0.15), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')
plt.gca().annotate(text_str_3, fontsize=12, xy=(0, -0.20), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')
plt.gca().annotate(text_str_2, fontsize=12, xy=(0.5, -0.15), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')
plt.gca().annotate(text_str_4, fontsize=12, xy=(0.5, -0.20), xycoords='axes fraction', xytext=(0, 0), textcoords='offset points', ha='left', va='top')

#Adjust layout to fit text
plt.subplots_adjust(bottom=0.2)

#save graph
plt.savefig(output_data + ".tiff", dpi=300)

#display graph
plt.show()

