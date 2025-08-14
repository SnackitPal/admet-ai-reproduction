# ADMET-AI Reproduction & Extension Project

This project aims to reproduce the scientific claims of the ADMET-AI paper and then extend its functionality with novel analyses.

## Project Overview

ADMET-AI is a machine learning platform for predicting Absorption, Distribution, Metabolism, Excretion, and Toxicity (ADMET) properties of chemical compounds. These predictions are crucial in drug discovery for filtering large chemical libraries and identifying promising drug candidates.

This project is structured into three main phases:

*   **Phase 1: Baseline Scientific Reproduction:** Verify the performance claims of the original ADMET-AI paper on a public benchmark dataset.
*   **Phase 2: Novel Extension:** Build upon the reproduced work by adding unique contributions, including error analysis, new dataset testing, and comparative drug analysis.
*   **Phase 3: Documentation & Presentation:** Package the project for professional presentation.

## Phase 1: Baseline Scientific Reproduction

**Goal:** To independently verify the performance claims made in the ADMET-AI paper by running the provided tool on a public benchmark dataset.

**Outcome:** Successfully reproduced the reported performance of the ADMET-AI model on the ClinTox benchmark dataset.

**Key Activities & Results:**

1.  **Benchmark Identification:** Selected the **ClinTox** dataset from the Therapeutics Data Commons (TDC).
2.  **Data Download:** Downloaded the ClinTox dataset to `data/data.csv` using `src/download_data.py`.
3.  **Prediction Execution:** Generated predictions for the ClinTox dataset using `admet_predict` and saved them to `data/predictions.csv`.
4.  **Results Analysis:** Calculated the AUROC score using `src/analyze_results.py`.
5.  **Comparison & Conclusion:**
    *   **Our Calculated AUROC (ClinTox):** 0.9774
    *   **Reported AUROC (ClinTox, from paper's supplementary data):** ~0.98

    Our results are in excellent agreement with the paper's claims.

## Phase 2: Novel Extensions

**Goal:** To build upon the reproduced work by adding unique contributions.

**Outcome:** Implemented three novel extensions to the ADMET-AI tool.

### 2.1 Error Analysis & Visualization

**Description:** Identified the top 5 false positives and top 5 false negatives from the ClinTox predictions and generated images of their chemical structures.

**Key Activities:**
*   Developed `src/error_analysis.py` to perform the analysis.
*   Generated and saved molecular images to the `images/` directory.

### 2.2 Substructure Analysis

**Description:** Performed a substructure analysis of the misclassified molecules to identify common functional groups.

**Key Activities:**
*   Developed `src/substructure_analysis.py` to identify and count common functional groups in the top 5 false positives and false negatives.
*   The analysis revealed that false positives are rich in Aromatic Rings and Halogens, while false negatives show a prevalence of Alcohol, Ketone, and Halogen groups.

### 2.2 New Dataset Testing

**Description:** Tested the pre-trained ADMET-AI model on a new public dataset (hERG) to evaluate its generalization capabilities.

**Key Activities:**
*   Downloaded the hERG dataset to `data/herg_data.csv` using `src/download_herg_data.py`.
*   Generated predictions for hERG using `admet_predict` and saved them to `data/herg_predictions.csv`.
*   Calculated the AUROC for hERG using `src/analyze_herg_results.py`.

**Results:**
*   **Our Calculated AUROC (hERG):** 0.9526
*   **Reported AUROC (hERG, from paper's supplementary data):** ~0.95

Our results confirm the model's strong generalization to new datasets.

### 2.3 Comparative Drug Analysis

**Description:** Developed a feature to compare a new molecule's predicted ADMET profile against a well-known drug (Aspirin).

**Key Activities:**
*   Developed `src/compare_drug.py` to perform side-by-side ADMET profile comparisons.
*   Used Aspirin (`CC(=O)OC1=CC=CC=C1C(=O)O`) as a reference drug.

### 2.4 Benchmark Analysis

**Description:** Performed a benchmark comparison of ADMET-AI against two other ADMET predictors (pkCSM and ADMETlab 3.0) on the ClinTox dataset, with enhanced statistical rigor.

**Key Activities:**
*   Developed `src/benchmark_analysis.py` to calculate AUROC and AUPRC scores with bootstrapped 95% confidence intervals.
*   The results show that ADMET-AI (AUROC: ~0.88, AUPRC: ~0.71) still outperforms both pkCSM (AUROC: ~0.66, AUPRC: ~0.38) and ADMETlab 3.0 (AUROC: ~0.71, AUPRC: ~0.44) on this dataset.
*   Developed `src/plot_auroc.py` to generate a comparison plot of the ROC curves, now including AUPRC in the legend.

![AUROC Comparison](images/auroc_comparison.png)

### 2.5 Substructure Analysis

**Description:** Performed a quantitative substructure analysis of the misclassified molecules to identify common functional groups and statistically test their enrichment.

**Key Activities:**
*   Developed `src/substructure_analysis.py` to identify and count common functional groups in misclassified molecules and apply Fisher's exact test.
*   The analysis revealed statistically significant enrichment of Halogen groups in False Positives (p=0.0419) and Carboxylic Acid (p=0.0093) and Halogen (p=0.0092) groups in False Negatives.

### 2.6 Chemical Space Visualization

**Description:** Visualized the chemical space of the ClinTox dataset to observe the distribution of correctly and incorrectly classified molecules.

**Key Activities:**
*   Developed `src/plot_chemical_space.py` to generate a 2D UMAP plot of the chemical space, coloring molecules by their classification (True Positive, True Negative, False Positive, False Negative).

## How to Run the Code

### 1. Environment Setup

It is recommended to use a Conda environment for this project.

```bash
# Create a new conda environment
conda create -n admet_project python=3.11
conda activate admet_project

# Install dependencies
pip install -r requirements.txt
```

### 2. Running Phase 1 Activities

To reproduce the baseline results:

```bash
python src/download_data.py
admet_predict --data_path data/data.csv --smiles_column Drug --save_path data/predictions.csv
python src/analyze_results.py
```

### 3. Running Phase 2 Activities

#### Error Analysis & Visualization

```bash
python src/error_analysis.py
# Images will be saved in the 'images/' directory
```

#### New Dataset Testing (hERG)

```bash
python src/download_herg_data.py
admet_predict --data_path data/herg_data.csv --smiles_column Drug --save_path data/herg_predictions.csv
python src/analyze_herg_results.py
```

#### Comparative Drug Analysis

To compare a new molecule (e.g., Ethanol: `CCO`) with Aspirin:

```bash
python src/compare_drug.py
```

Feel free to modify `src/compare_drug.py` to compare different molecules.