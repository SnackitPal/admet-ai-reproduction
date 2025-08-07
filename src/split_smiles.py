import os

def split_smiles_file(input_path='data/clintox_smiles.txt', output_dir='data/smiles_chunks', chunk_size=100):
    """
    Splits a file of SMILES strings into smaller chunk files.

    Args:
        input_path (str): Path to the input file with one SMILES per line.
        output_dir (str): Directory to save the chunk files.
        chunk_size (int): Number of SMILES strings per chunk file.
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    with open(input_path, 'r') as f:
        lines = f.readlines()

    for i in range(0, len(lines), chunk_size):
        chunk = lines[i:i + chunk_size]
        output_filename = os.path.join(output_dir, f"clintox_smiles_part_{(i // chunk_size) + 1}.txt")
        with open(output_filename, 'w') as chunk_file:
            chunk_file.writelines(chunk)
        print(f"Created {output_filename} with {len(chunk)} SMILES.")

if __name__ == "__main__":
    split_smiles_file()
