
import os
import pandas as pd
from admet_ai.admet_predict import admet_predict
from pathlib import Path

def compare_drug_profiles(new_molecule_smiles: str, reference_drug_smiles: str, reference_drug_name: str = "Aspirin"):
    """
    Compares the ADMET profiles of a new molecule and a reference drug.
    """
    # Create a DataFrame for the new molecule
    new_mol_df = pd.DataFrame({'Drug': [new_molecule_smiles]})

    # Create a DataFrame for the reference drug
    ref_drug_df = pd.DataFrame({'Drug': [reference_drug_smiles]})

    # Save DataFrames to temporary CSV files
    new_mol_temp_path = 'data/new_molecule_temp.csv'
    ref_drug_temp_path = 'data/reference_drug_temp.csv'
    new_mol_df.to_csv(new_mol_temp_path, index=False)
    ref_drug_df.to_csv(ref_drug_temp_path, index=False)

    # Get predictions for the new molecule
    admet_predict(data_path=new_mol_temp_path, smiles_column='Drug', save_path=Path('data/new_molecule_predictions.csv'))
    new_mol_predictions = pd.read_csv('data/new_molecule_predictions.csv')

    # Get predictions for the reference drug
    admet_predict(data_path=ref_drug_temp_path, smiles_column='Drug', save_path=Path('data/reference_drug_predictions.csv'))
    ref_drug_predictions = pd.read_csv('data/reference_drug_predictions.csv')

    # Clean up temporary files
    os.remove(new_mol_temp_path)
    os.remove(ref_drug_temp_path)
    os.remove('data/new_molecule_predictions.csv')
    os.remove('data/reference_drug_predictions.csv')

    print(f"\nADMET Profile Comparison: New Molecule vs. {reference_drug_name}")
    print("------------------------------------------------------------------")

    # Display predictions side-by-side
    comparison_df = pd.concat([new_mol_predictions.iloc[0], ref_drug_predictions.iloc[0]], axis=1)
    comparison_df.columns = ["New Molecule", reference_drug_name]
    print(comparison_df.to_string())

if __name__ == "__main__":
    # Example usage: Replace with your new molecule's SMILES
    new_smiles = "CCO"
    aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"

    compare_drug_profiles(new_smiles, aspirin_smiles, "Aspirin")
