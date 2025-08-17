
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import os

def error_analysis():
    """
    Performs error analysis by identifying the top 5 false positives and false negatives,
    and generates images of their chemical structures.
    """
    # Load the original data
    true_data = pd.read_csv('admet_predict/data/data.csv')

    # Load the predictions
    predictions = pd.read_csv('admet_predict/data/predictions.csv')

    # Merge the dataframes on the 'Drug' column
    df = pd.merge(true_data, predictions, on='Drug', how='inner')

    print("\nMerged DataFrame columns:", df.columns)
    print("\nMerged DataFrame head:\n", df.head())

    # Ensure 'Y' and 'ClinTox' columns are numeric
    # Ensure 'Y_x' (true labels) and 'ClinTox' (predicted probabilities) columns are numeric
    df['Y_x'] = pd.to_numeric(df['Y_x'])
    df['ClinTox'] = pd.to_numeric(df['ClinTox'])

    # Calculate the error (difference between predicted probability and true label)
    df['error'] = abs(df['Y_x'] - df['ClinTox'])

    # Separate false positives and false negatives
    false_positives = df[(df['Y_x'] == 0) & (df['ClinTox'] > 0.5)].nlargest(5, 'error')
    false_negatives = df[(df['Y_x'] == 1) & (df['ClinTox'] < 0.5)].nlargest(5, 'error')

    print("Top 5 False Positives:")
    print(false_positives[['Drug', 'Y_x', 'ClinTox', 'error']])

    print("\nTop 5 False Negatives:")
    print(false_negatives[['Drug', 'Y_x', 'ClinTox', 'error']])

    # Create the images directory if it doesn't exist
    if not os.path.exists('admet_predict/images'):
        os.makedirs('admet_predict/images')

    # Generate and save images
    for index, row in false_positives.iterrows():
        mol = Chem.MolFromSmiles(row['Drug'])
        Draw.MolToFile(mol, f"admet_predict/images/fp_{index}.png")

    for index, row in false_negatives.iterrows():
        mol = Chem.MolFromSmiles(row['Drug'])
        Draw.MolToFile(mol, f"admet_predict/images/fn_{index}.png")

    print("\nImages of the top 5 false positives and false negatives have been saved to the 'admet_predict/images/' directory.")

if __name__ == "__main__":
    error_analysis()
