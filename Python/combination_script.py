import pandas as pd
import os

# Prompt user for directory containing TSV files
tsv_directory = input("Enter the path to the directory containing TSV files: ").strip()

# Check if the directory exists
if not os.path.isdir(tsv_directory):
    print(f"The directory '{tsv_directory}' does not exist.")
    exit()

# Prompt for the output file name
output_workbook = input("Enter the output Excel file name (e.g., combined_workbook.xlsx): ").strip()

# Add .xlsx extension if not provided
if not output_workbook.endswith(".xlsx"):
    output_workbook += ".xlsx"

# Initialize an Excel writer
with pd.ExcelWriter(output_workbook, engine='openpyxl') as writer:
    # Process TSV files in the directory
    for file in os.listdir(tsv_directory):
        if file.endswith(".tsv"):  # Only process .tsv files
            file_path = os.path.join(tsv_directory, file)
            sheet_name = os.path.splitext(file)[0]  # Use the file name (without extension) as the sheet name
            df = pd.read_csv(file_path, sep="\t")
            df.to_excel(writer, index=False, sheet_name=sheet_name[:31])  # Limit sheet name to 31 characters

print(f"Workbook saved as '{output_workbook}'")