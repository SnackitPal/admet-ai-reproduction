import pandas as pd
from rdkit import Chem
from collections import defaultdict
from scipy.stats import fisher_exact

def analyze_substructures_presence(smiles_list):
    """
    Analyzes a list of SMILES strings to determine the presence of common functional groups.

    Args:
        smiles_list (list): A list of SMILES strings.

    Returns:
        dict: A dictionary where keys are functional group names and values are lists of booleans
              (True if present, False if not) for each molecule in the input smiles_list.
    """
    functional_groups = {
        "Alcohol": "[#6][OD1H]",
        "Amine": "[NX3;H2;!$(NC=O)]",
        "Aromatic Ring": "a",
        "Carboxylic Acid": "[CX3](=O)[OX2H1]",
        "Ester": "[#6]C(=O)O[#6]",
        "Ether": "[OD2]([#6])[#6]",
        "Ketone": "[#6](=O)[#6]",
        "Nitro": "[$([NX3](=O)=O),-]([#6])",
        "Sulfonamide": "S(=O)(=O)N",
        "Halogen (F, Cl, Br, I)": "[F,Cl,Br,I]",
    }

    presence_data = defaultdict(list)
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # If SMILES is invalid, assume no functional groups are present
            for name in functional_groups:
                presence_data[name].append(False)
            continue

        for name, smarts in functional_groups.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                presence_data[name].append(False)
                continue
            presence_data[name].append(mol.HasSubstructMatch(pattern))
    
    return presence_data

# Load the true labels and predicted probabilities
true_labels_df = pd.read_csv('admet_predict/data/data.csv', header=0)
predictions_df = pd.read_csv('admet_predict/data/predictions.csv', header=0)

# Merge the dataframes on the 'Drug' column
merged_df = pd.merge(true_labels_df, predictions_df, on='Drug', how='inner')

# Ensure 'Y_x' and 'ClinTox' are numeric
merged_df['Y_x'] = pd.to_numeric(merged_df['Y_x'], errors='coerce')
merged_df['ClinTox'] = pd.to_numeric(merged_df['ClinTox'], errors='coerce')

# Drop rows with NaN values that might result from coercion
merged_df.dropna(subset=['Y_x', 'ClinTox'], inplace=True)

# --- Define classification groups ---
# True Positives: Y=1, Pred > 0.5
tp_df = merged_df[(merged_df['Y_x'] == 1) & (merged_df['ClinTox'] > 0.5)]
tp_smiles = tp_df['Drug'].tolist()

# True Negatives: Y=0, Pred < 0.5
tn_df = merged_df[(merged_df['Y_x'] == 0) & (merged_df['ClinTox'] < 0.5)]
tn_smiles = tn_df['Drug'].tolist()

# False Positives: Y=0, Pred > 0.5
fp_df = merged_df[(merged_df['Y_x'] == 0) & (merged_df['ClinTox'] > 0.5)]
fp_smiles = fp_df['Drug'].tolist()

# False Negatives: Y=1, Pred < 0.5
fn_df = merged_df[(merged_df['Y_x'] == 1) & (merged_df['ClinTox'] < 0.5)]
fn_smiles = fn_df['Drug'].tolist()

# --- Perform substructure analysis for all groups ---
tp_presence = analyze_substructures_presence(tp_smiles)
tn_presence = analyze_substructures_presence(tn_smiles)
fp_presence = analyze_substructures_presence(fp_smiles)
fn_presence = analyze_substructures_presence(fn_smiles)

# --- Perform Fisher's Exact Test and print results ---
print("---" + "-" * 40)
print("--- Statistical Analysis of Functional Group Enrichment ---")
print("\nComparing False Positives (FP) vs. True Negatives (TN):")
print("{:<20} | {:<10} | {:<10}".format('Functional Group', 'FP Count', 'TN Count'))
print("-" * 45)
for fg_name in tp_presence.keys(): # Use keys from any presence dict
    fp_has_fg = sum(fp_presence[fg_name])
    fp_no_fg = len(fp_presence[fg_name]) - fp_has_fg
    tn_has_fg = sum(tn_presence[fg_name])
    tn_no_fg = len(tn_presence[fg_name]) - tn_has_fg

    table = [[fp_has_fg, fp_no_fg],
             [tn_has_fg, tn_no_fg]]
    
    if all(sum(row) > 0 for row in table) and all(sum(col) > 0 for col in zip(*table)):
        oddsratio, pvalue = fisher_exact(table)
        print(f"{fg_name:<20} | {fp_has_fg:<10} | {tn_has_fg:<10} | p-value: {pvalue:.4f}")
    else:
        print(f"{fg_name:<20} | {fp_has_fg:<10} | {tn_has_fg:<10} | Not enough data for test")

print("\nComparing False Negatives (FN) vs. True Positives (TP):")
print("{:<20} | {:<10} | {:<10}".format('Functional Group', 'FN Count', 'TP Count'))
print("-" * 45)
for fg_name in tp_presence.keys():
    fn_has_fg = sum(fn_presence[fg_name])
    fn_no_fg = len(fn_presence[fg_name]) - fn_has_fg
    tp_has_fg = sum(tp_presence[fg_name])
    tp_no_fg = len(tp_presence[fg_name]) - tp_has_fg

    table = [[fn_has_fg, fn_no_fg],
             [tp_has_fg, tp_no_fg]]

    if all(sum(row) > 0 for row in table) and all(sum(col) > 0 for col in zip(*table)):
        oddsratio, pvalue = fisher_exact(table)
        print(f"{fg_name:<20} | {fn_has_fg:<10} | {tp_has_fg:<10} | p-value: {pvalue:.4f}")
    else:
        print(f"{fg_name:<20} | {fn_has_fg:<10} | {tp_has_fg:<10} | Not enough data for test")
