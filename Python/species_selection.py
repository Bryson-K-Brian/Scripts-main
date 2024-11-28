import pandas as pd

# Input Excel file
input_workbook = input("Enter the path to the Excel file: ").strip()

# Column name to filter
target_column = input("Enter the column name to filter: ").strip()

# Filter condition (e.g., keep rows where the column value equals a specific value)
filter_value = input("Enter the value to filter by: ").strip()

# Output file name
output_workbook = input("Enter the output Excel file name (e.g., filtered_workbook.xlsx): ").strip()

# Add .xlsx extension if not provided
if not output_workbook.endswith(".xlsx"):
    output_workbook += ".xlsx"

# Open the workbook
with pd.ExcelWriter(output_workbook, engine='openpyxl') as writer:
    # Read all worksheets into a dictionary of DataFrames
    xls = pd.ExcelFile(input_workbook)
    for sheet_name in xls.sheet_names:
        # Read the sheet into a DataFrame
        df = pd.read_excel(input_workbook, sheet_name=sheet_name)
        
        # Check if the target column exists in this sheet
        if target_column in df.columns:
            # Apply the filter (keep rows where the column equals the filter value)
            filtered_df = df[df[target_column] == filter_value]
        else:
            # If the column doesn't exist, keep the sheet as is
            filtered_df = df

        # Write the filtered DataFrame to the output workbook
        filtered_df.to_excel(writer, index=False, sheet_name=sheet_name)

print(f"Filtered workbook saved as '{output_workbook}'")