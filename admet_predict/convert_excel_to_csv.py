import pandas as pd
import os

def convert_excel_sheets_to_csv(excel_file_path):
    """
    Reads an Excel file, converts each sheet to a separate CSV file,
    and saves them in a new subdirectory.

    Args:
        excel_file_path (str): The absolute path to the Excel file.
    """
    if not os.path.exists(excel_file_path):
        print(f"Error: Excel file not found at '{excel_file_path}'")
        return

    # Get the directory of the Excel file
    excel_dir = os.path.dirname(excel_file_path)
    
    # Create the output directory
    output_dir = os.path.join(excel_dir, "excel_sheets_csv")
    os.makedirs(output_dir, exist_ok=True)

    print(f"Converting Excel file: '{excel_file_path}'")
    print(f"Saving CSVs to: '{output_dir}'")

    try:
        # Read all sheets from the Excel file
        xls = pd.ExcelFile(excel_file_path)
        
        for sheet_name in xls.sheet_names:
            df = xls.parse(sheet_name)
            
            # Sanitize sheet name for filename
            sanitized_sheet_name = "".join([c for c in sheet_name if c.isalnum() or c in [' ', '_', '-']]).replace(' ', '_')
            
            # Construct output filename
            output_filename = f"Supplementary_Table_1_{sanitized_sheet_name}.csv"
            output_file_path = os.path.join(output_dir, output_filename)
            
            df.to_csv(output_file_path, index=False)
            print(f"  - Converted sheet '{sheet_name}' to '{output_filename}'")
            
    except Exception as e:
        print(f"An error occurred during conversion: {e}")

# --- Script execution ---
excel_file = "/home/sanket/admet-ai-reproduction/btae416_supplementary_data/Supplementary Table 1.xlsx"
convert_excel_sheets_to_csv(excel_file)
