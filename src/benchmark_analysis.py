import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score, matthews_corrcoef
import numpy as np
from scipy.stats import wilcoxon

def load_and_prepare_data():
    """Loads, prepares, and merges the data from the different sources."""
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

    # Ensure prediction columns are numeric
    for model_name in ['ADMET-AI', 'pkCSM', 'ADMETlab']:
        final_df[model_name] = pd.to_numeric(final_df[model_name], errors='coerce')

    # Drop any rows where data might be missing to be safe
    final_df.dropna(inplace=True)

    return final_df

def bootstrap_metrics(y_true, y_pred, n_bootstraps=1000):
    """Calculates bootstrapped metrics (AUROC, AUPRC, MCC)."""
    auroc_scores = []
    auprc_scores = []
    mcc_scores = []
    n_samples = len(y_true)
    for _ in range(n_bootstraps):
        indices = np.random.choice(n_samples, n_samples, replace=True)
        y_true_sample = y_true.iloc[indices]
        y_pred_sample = y_pred.iloc[indices]
        
        if len(np.unique(y_true_sample)) < 2:
            continue

        auroc_scores.append(roc_auc_score(y_true_sample, y_pred_sample))
        auprc_scores.append(average_precision_score(y_true_sample, y_pred_sample))
        
        # For MCC, we need binary predictions
        y_pred_binary = (y_pred_sample > 0.5).astype(int)
        mcc_scores.append(matthews_corrcoef(y_true_sample, y_pred_binary))

    return auroc_scores, auprc_scores, mcc_scores

def run_bootstrap_analysis(df, models, n_bootstraps=1000):
    """Runs the bootstrap analysis for all models."""
    results = {}
    bootstrap_scores = {}
    for model_name in models:
        auroc_scores, auprc_scores, mcc_scores = bootstrap_metrics(df['Y'], df[model_name], n_bootstraps)
        
        results[model_name] = {
            'AUROC_mean': np.mean(auroc_scores),
            'AUROC_ci_lower': np.percentile(auroc_scores, 2.5),
            'AUROC_ci_upper': np.percentile(auroc_scores, 97.5),
            'AUPRC_mean': np.mean(auprc_scores),
            'AUPRC_ci_lower': np.percentile(auprc_scores, 2.5),
            'AUPRC_ci_upper': np.percentile(auprc_scores, 97.5),
            'MCC_mean': np.mean(mcc_scores),
            'MCC_ci_lower': np.percentile(mcc_scores, 2.5),
            'MCC_ci_upper': np.percentile(mcc_scores, 97.5),
        }
        bootstrap_scores[model_name] = auroc_scores
    return results, bootstrap_scores

def calculate_statistical_significance(bootstrap_scores, models):
    """Calculates statistical significance between models."""
    p_values = {}
    for i in range(len(models)):
        for j in range(i + 1, len(models)):
            model1 = models[i]
            model2 = models[j]
            stat, p_value = wilcoxon(bootstrap_scores[model1], bootstrap_scores[model2])
            p_values[f'{model1}_vs_{model2}'] = p_value
    return p_values

def print_and_save_results(results, p_values):
    """Prints and saves the results table."""
    print("--- ClinTox Benchmark Comparison with Bootstrapped CIs ---")
    print("\n{:<15} | {:<10} | {:<10} | {:<10} | {:<10} | {:<10} | {:<10}".format('Model', 'AUROC', 'AUROC CI', 'AUPRC', 'AUPRC CI', 'MCC', 'MCC CI'))
    print("-" * 100)
    
    results_df = pd.DataFrame(columns=['Model', 'AUROC', 'AUROC CI', 'AUPRC', 'AUPRC CI', 'MCC', 'MCC CI'])

    for model, metrics in results.items():
        auroc_ci = f"({metrics['AUROC_ci_lower']:.2f}-{metrics['AUROC_ci_upper']:.2f})"
        auprc_ci = f"({metrics['AUPRC_ci_lower']:.2f}-{metrics['AUPRC_ci_upper']:.2f})"
        mcc_ci = f"({metrics['MCC_ci_lower']:.2f}-{metrics['MCC_ci_upper']:.2f})"
        
        print("{:<15} | {:<10.4f} | {:<10} | {:<10.4f} | {:<10} | {:<10.4f} | {:<10}".format(
            model,
            metrics['AUROC_mean'],
            auroc_ci,
            metrics['AUPRC_mean'],
            auprc_ci,
            metrics['MCC_mean'],
            mcc_ci
        ))
        
        results_df.loc[len(results_df)] = [model, metrics['AUROC_mean'], auroc_ci, metrics['AUPRC_mean'], auprc_ci, metrics['MCC_mean'], mcc_ci]

    print("\n--- Statistical Significance (p-values) ---")
    for pair, p_value in p_values.items():
        print(f"{pair}: {p_value:.4f}")

    print("\n----------------------------------------------------------------------")
    
    results_df.to_csv('data/benchmark_results.csv', index=False)
    print("\nResults saved to data/benchmark_results.csv")


def main():
    """Main function to run the benchmark analysis."""
    models_to_compare = ['ADMET-AI', 'pkCSM', 'ADMETlab']
    
    final_df = load_and_prepare_data()
    results, bootstrap_scores = run_bootstrap_analysis(final_df, models_to_compare)
    p_values = calculate_statistical_significance(bootstrap_scores, models_to_compare)
    print_and_save_results(results, p_values)

if __name__ == "__main__":
    main()
