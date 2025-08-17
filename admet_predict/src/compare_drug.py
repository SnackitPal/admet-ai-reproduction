import pandas as pd
from admet_ai import ADMETModel

def compare_admet_profiles(new_molecule_smiles, reference_drug_smiles, new_molecule_name="New Molecule", reference_drug_name="Aspirin"):
    """
    Generates and compares the ADMET profiles of two molecules.

    Args:
        new_molecule_smiles (str): SMILES string of the new molecule.
        reference_drug_smiles (str): SMILES string of the reference drug.
        new_molecule_name (str, optional): Name of the new molecule. Defaults to "New Molecule".
        reference_drug_name (str, optional): Name of the reference drug. Defaults to "Aspirin".
    """
    # Initialize the ADMET-AI model
    model = ADMETModel()

    # Create a list of SMILES strings to predict
    smiles_list = [new_molecule_smiles, reference_drug_smiles]

    # Generate predictions
    preds_df = model.predict(smiles=smiles_list)

    # Create a mapping from SMILES to desired names
    smiles_to_name = {
        new_molecule_smiles: new_molecule_name,
        reference_drug_smiles: reference_drug_name
    }

    # Rename the index using the mapping
    preds_df.rename(index=smiles_to_name, inplace=True)

    # Transpose the dataframe for side-by-side comparison
    comparison_df = preds_df.transpose()

    # Print the comparison
    print(f"ADMET Profile Comparison: {new_molecule_name} vs. {reference_drug_name}")
    print("-" * (len(new_molecule_name) + len(reference_drug_name) + 30))
    print(comparison_df[[new_molecule_name, reference_drug_name]])


if __name__ == "__main__":
    # --- Molecule Definitions ---
    # You can modify these SMILES strings to compare different molecules
    new_molecule_smiles = "CCO"  # Ethanol
    reference_drug_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin

    compare_admet_profiles(new_molecule_smiles, reference_drug_smiles)