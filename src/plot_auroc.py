import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, average_precision_score

# --- Load all necessary files ---
true_labels_df = pd.read_csv('data/data.csv')
admet_ai_preds = pd.read_csv('data/predictions.csv')
pkcsm_preds = pd.read_csv('data/pkcsm_predictions_combined.csv')
admetlab_preds = pd.read_csv('data/admetlab_predictions.csv')

# --- Prepare the data ---
# 1. True Labels
true_labels = true_labels_df['Y']

# 2. ADMET-AI Predictions
admet_ai_scores = admet_ai_preds['ClinTox']

# 3. pkCSM Predictions
pkcsm_preds['pkCSM_prob'] = pkcsm_preds['Hepatotoxicity'].apply(lambda x: 1 if x == 'Yes' else 0)
pkcsm_scores = pkcsm_preds['pkCSM_prob']

# 4. ADMETlab 2.0 Predictions
admetlab_scores = admetlab_preds['H-HT']


# --- Align data ---
final_df = pd.DataFrame({
    'Y': true_labels,
    'ADMET-AI': admet_ai_scores,
    'pkCSM': pkcsm_scores,
    'ADMETlab': admetlab_scores
}).dropna()


# --- Calculate ROC curve and AUC for each model ---
plt.figure(figsize=(8, 6))

models = {
    'ADMET-AI': final_df['ADMET-AI'],
    'pkCSM': final_df['pkCSM'],
    'ADMETlab 2.0': final_df['ADMETlab']
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
plt.savefig('images/auroc_comparison.png', dpi=300)

print("AUROC comparison plot saved to images/auroc_comparison.png")