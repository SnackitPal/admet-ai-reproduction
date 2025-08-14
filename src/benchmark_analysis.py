import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score
import numpy as np

# --- Load all necessary files ---
true_labels_df = pd.read_csv('data/data.csv')
admet_ai_preds = pd.read_csv('data/predictions.csv')
pkcsm_preds = pd.read_csv('data/pkcsm_predictions_combined.csv')
admetlab_preds = pd.read_csv('data/admetlab_predictions_cleaned.csv')

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

# Drop any rows where data might be missing to be safe
final_df.dropna(inplace=True)

# --- Bootstrapping function ---
def bootstrap_metrics(y_true, y_pred, n_bootstraps=1000):
    auroc_scores = []
    auprc_scores = []
    n_samples = len(y_true)
    for _ in range(n_bootstraps):
        indices = np.random.choice(n_samples, n_samples, replace=True)
        y_true_sample = y_true.iloc[indices]
        y_pred_sample = y_pred.iloc[indices]
        
        # Ensure there's at least one positive and one negative sample in the bootstrap sample
        if len(np.unique(y_true_sample)) < 2:
            continue

        auroc_scores.append(roc_auc_score(y_true_sample, y_pred_sample))
        auprc_scores.append(average_precision_score(y_true_sample, y_pred_sample))
    return auroc_scores, auprc_scores

# --- Calculate AUROC and AUPRC with bootstrapping for each model ---
results = {}
for model_name in ['ADMET-AI', 'pkCSM', 'ADMETlab']:
    # Ensure the prediction column is numeric
    final_df[model_name] = pd.to_numeric(final_df[model_name], errors='coerce')
    # Drop rows where coercion might have created NaNs
    temp_df = final_df[['Y', model_name]].dropna()
    
    auroc_scores, auprc_scores = bootstrap_metrics(temp_df['Y'], temp_df[model_name])
    
    results[model_name] = {
        'AUROC_mean': np.mean(auroc_scores),
        'AUROC_ci_lower': np.percentile(auroc_scores, 2.5),
        'AUROC_ci_upper': np.percentile(auroc_scores, 97.5),
        'AUPRC_mean': np.mean(auprc_scores),
        'AUPRC_ci_lower': np.percentile(auprc_scores, 2.5),
        'AUPRC_ci_upper': np.percentile(auprc_scores, 97.5),
    }

# --- Print the results table ---
print("--- ClinTox Benchmark Comparison with Bootstrapped CIs ---")
print("\n{:<15} | {:<10} | {:<10} | {:<10} | {:<10}".format('Model', 'AUROC', 'AUROC CI', 'AUPRC', 'AUPRC CI'))
print("-" * 70)
for model, metrics in results.items():
    print("{:<15} | {:<10.4f} | ({:<.2f}-{:<.2f}) | {:<10.4f} | ({:<.2f}-{:<.2f})".format(
        model,
        metrics['AUROC_mean'],
        metrics['AUROC_ci_lower'],
        metrics['AUROC_ci_upper'],
        metrics['AUPRC_mean'],
        metrics['AUPRC_ci_lower'],
        metrics['AUPRC_ci_upper']
    ))
print("\n----------------------------------------------------------------------")