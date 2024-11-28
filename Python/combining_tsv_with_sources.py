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

# List to hold all the DataFrames
combined_df = []

# Process TSV files in the directory
for file in os.listdir(tsv_directory):
    if file.endswith(".tsv"):  # Only process .tsv files
        file_path = os.path.join(tsv_directory, file)
        print(f"Processing file: {file}")
        
        try:
            # Read the TSV file into a DataFrame
            df = pd.read_csv(file_path, sep="\t", skip_blank_lines=True)
            
            # Drop rows where all columns are NaN
            df.dropna(how="all", inplace=True)
            
            # Ensure unique column names by appending a suffix to duplicates
            df.columns = pd.Index([f"{col}.{i}" if df.columns.tolist().count(col) > 1 else col 
                                   for i, col in enumerate(df.columns)])
            
            # Add a new column with the sheet name (or file name without extension)
            df['Source_Sheet'] = os.path.splitext(file)[0]
            
            # Reset index to avoid duplicate indices
            df.reset_index(drop=True, inplace=True)
            
            # Append the DataFrame to the list
            combined_df.append(df)
        except Exception as e:
            print(f"Error reading '{file_path}': {e}")
            continue

# Concatenate all DataFrames into a single one
if combined_df:
    try:
        final_df = pd.concat(combined_df, ignore_index=True)
        # Save the combined DataFrame into the output Excel file
        final_df.to_excel(output_workbook, index=False)
        print(f"Workbook saved as '{output_workbook}'")
    except Exception as e:
        print(f"Error during concatenation or saving: {e}")
else:
    print("No data was combined. Check your files for inconsistencies.")