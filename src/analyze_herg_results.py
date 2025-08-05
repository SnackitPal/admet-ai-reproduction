
import pandas as pd
from sklearn.metrics import roc_auc_score

def analyze_herg_results():
    """
    Analyzes the hERG prediction results by calculating the AUROC score.
    """
    # Load the original hERG data to get the true labels
    true_data = pd.read_csv('data/herg_data.csv')
    true_labels = true_data['Y']

    # Load the hERG predictions
    predictions = pd.read_csv('data/herg_predictions.csv')
    # The relevant column for hERG is 'hERG'
    predicted_probabilities = predictions['hERG']

    # Calculate the AUROC score
    auroc = roc_auc_score(true_labels, predicted_probabilities)

    print(f"AUROC Score for hERG: {auroc:.4f}")

if __name__ == "__main__":
    analyze_herg_results()
