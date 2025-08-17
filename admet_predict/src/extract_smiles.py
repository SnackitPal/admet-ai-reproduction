import pandas as pd

# Load the ClinTox dataset
df = pd.read_csv('admet_predict/data/data.csv')

# Extract the SMILES strings from the 'Drug' column
smiles_list = df['Drug'].tolist()

# Write the SMILES strings to a text file
with open('admet_predict/data/clintox_smiles.txt', 'w') as f:
    for smiles in smiles_list:
        f.write(f"{smiles}\n")

print("SMILES strings extracted and saved to admet_predict/data/clintox_smiles.txt")
