import os
import pandas as pd

# Prompt user for directory containing TSV files
tsv_directory = input("Enter the path to the directory containing TSV files: ").strip()

# Check if the directory exists
if not os.path.isdir(tsv_directory):
    print(f"The directory '{tsv_directory}' does not exist.")
    exit()

# Define the new header names
headers = ["Percentage", "Rollover_counts", "Exact_Counts", "Taxonomy_level", "Tax_ID", "Taxa_Name"]  # Replace with your desired column names

# Process each TSV file in the directory
for file in os.listdir(tsv_directory):
    if file.endswith(".tsv"):
        file_path = os.path.join(tsv_directory, file)
        print(f"Processing file: {file}")
        
        try:
            # Read the TSV file into a DataFrame without headers
            df = pd.read_csv(file_path, sep="\t", header=None)  # Treat file as if it has no headers
            
            # Assign new header names to the DataFrame
            df.columns = headers[:len(df.columns)]  # Ensure we don't exceed the number of actual columns
            
            # Save the updated DataFrame back to the TSV file
            df.to_csv(file_path, sep="\t", index=False)
            print(f"Updated headers for file: {file}")
        except Exception as e:
            print(f"Error processing '{file_path}': {e}")