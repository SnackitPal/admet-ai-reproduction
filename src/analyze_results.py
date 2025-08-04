
import pandas as pd
from sklearn.metrics import roc_auc_score

def analyze_results():
    """
    Analyzes the prediction results by calculating the AUROC score.
    """
    # Load the original data to get the true labels
    true_data = pd.read_csv('data/data.csv')
    true_labels = true_data['Y']

    # Load the predictions
    predictions = pd.read_csv('data/predictions.csv')
    # The relevant column for ClinTox is 'ClinTox'
    predicted_probabilities = predictions['ClinTox']

    # Calculate the AUROC score
    auroc = roc_auc_score(true_labels, predicted_probabilities)

    print(f"AUROC Score for ClinTox: {auroc:.4f}")

if __name__ == "__main__":
    analyze_results()
