import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, average_precision_score

# --- Load all necessary files ---
true_labels_df = pd.read_csv('admet_predict/data/data.csv')
admet_ai_preds = pd.read_csv('admet_predict/data/predictions.csv')
pkcsm_preds = pd.read_csv('admet_predict/data/pkcsm_predictions_combined.csv')
admetlab_preds = pd.read_csv('admet_predict/data/admetlab_predictions_cleaned.csv')

# --- Prepare the data ---
# Rename SMILES columns for consistency
true_labels_df.rename(columns={'Drug': 'SMILES'}, inplace=True)
admet_ai_preds.rename(columns={'Drug': 'SMILES'}, inplace=True)
admetlab_preds.rename(columns={'smiles': 'SMILES'}, inplace=True)

# Select relevant columns
true_labels = true_labels_df[['SMILES', 'Y']]
admet_ai_values = admet_ai_preds[['SMILES', 'ClinTox']]
pkcsm_preds['pkCSM'] = pkcsm_preds['Hepatotoxicity'].apply(lambda x: 1 if x == 'Yes' else 0)
pkcsm_values = pkcsm_preds[['SMILES', 'pkCSM']]
admetlab_values = admetlab_preds[['SMILES', 'H-HT']]

# Merge dataframes on SMILES string
final_df = pd.merge(true_labels, admet_ai_values, on='SMILES', how='inner')
final_df = pd.merge(final_df, pkcsm_values, on='SMILES', how='inner')
final_df = pd.merge(final_df, admetlab_values, on='SMILES', how='inner')

# Rename columns for the analysis loop
final_df.rename(columns={'ClinTox': 'ADMET-AI', 'H-HT': 'ADMETlab'}, inplace=True)

# Ensure prediction columns are numeric
final_df['ADMET-AI'] = pd.to_numeric(final_df['ADMET-AI'], errors='coerce')
final_df['pkCSM'] = pd.to_numeric(final_df['pkCSM'], errors='coerce')
final_df['ADMETlab'] = pd.to_numeric(final_df['ADMETlab'], errors='coerce')

final_df.dropna(inplace=True)


# --- Calculate ROC curve and AUC for each model ---
plt.figure(figsize=(8, 6))

models = {
    'ADMET-AI': final_df['ADMET-AI'],
    'pkCSM': final_df['pkCSM'],
    'ADMETlab 3.0': final_df['ADMETlab']
}

for name, scores in models.items():
    fpr, tpr, _ = roc_curve(final_df['Y'], scores)
    roc_auc = auc(fpr, tpr)
    auprc = average_precision_score(final_df['Y'], scores)
    plt.plot(fpr, tpr, lw=2, label=f'{name} (AUROC = {roc_auc:.2f}, AUPRC = {auprc:.2f})')

# --- Plotting ---
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve Comparison')
plt.legend(loc="lower right")
plt.grid(True)
plt.savefig('admet_predict/images/auroc_comparison.png', dpi=300)

print("AUROC comparison plot saved to admet_predict/images/auroc_comparison.png")