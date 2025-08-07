import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import umap
import matplotlib.pyplot as plt
import numpy as np

# --- Load data ---
true_labels_df = pd.read_csv('data/data.csv')
predictions_df = pd.read_csv('data/predictions.csv')

# Merge dataframes
merged_df = pd.merge(true_labels_df, predictions_df, on='Drug', how='inner')
merged_df['Y_x'] = pd.to_numeric(merged_df['Y_x'], errors='coerce')
merged_df['ClinTox'] = pd.to_numeric(merged_df['ClinTox'], errors='coerce')
merged_df.dropna(subset=['Y_x', 'ClinTox'], inplace=True)

# --- Classify molecules ---
def classify_molecule(row):
    true_label = row['Y_x']
    predicted_prob = row['ClinTox']
    
    if true_label == 1 and predicted_prob >= 0.5:
        return 'True Positive'
    elif true_label == 0 and predicted_prob < 0.5:
        return 'True Negative'
    elif true_label == 0 and predicted_prob >= 0.5:
        return 'False Positive'
    elif true_label == 1 and predicted_prob < 0.5:
        return 'False Negative'
    return 'Unknown'

merged_df['Classification'] = merged_df.apply(classify_molecule, axis=1)

# --- Generate Morgan Fingerprints ---
mols = [Chem.MolFromSmiles(s) for s in merged_df['Drug']]
valid_mols = [m for m in mols if m is not None]
valid_indices = [i for i, m in enumerate(mols) if m is not None]

fingerprints = [AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in valid_mols]
fingerprints_np = np.array(fingerprints)

# Filter merged_df to include only valid molecules
filtered_df = merged_df.iloc[valid_indices].copy()

# --- UMAP Dimensionality Reduction ---
reducer = umap.UMAP(n_components=2, random_state=42)
umap_embedding = reducer.fit_transform(fingerprints_np)

filtered_df['UMAP_x'] = umap_embedding[:, 0]
filtered_df['UMAP_y'] = umap_embedding[:, 1]

# --- Plotting ---
plt.figure(figsize=(10, 8))

colors = {
    'True Positive': '#2ca02c',  # Green
    'True Negative': '#1f77b4',  # Blue
    'False Positive': '#d62728', # Red
    'False Negative': '#ff7f0e'  # Orange
}

for classification, color in colors.items():
    subset = filtered_df[filtered_df['Classification'] == classification]
    plt.scatter(subset['UMAP_x'], subset['UMAP_y'], c=color, label=classification, alpha=0.7, s=10)

plt.title('Chemical Space Visualization of ClinTox Dataset (UMAP)')
plt.xlabel('UMAP Dimension 1')
plt.ylabel('UMAP Dimension 2')
plt.legend()
plt.grid(True)
plt.savefig('images/chemical_space_plot.png', dpi=300)

print("Chemical space plot saved to images/chemical_space_plot.png")
