# Project Progress Log

This document provides a detailed log of the ADMET-AI reproduction and extension project. Its purpose is to track our progress and serve as a persistent source of context.

---

## Phase 1: Baseline Scientific Reproduction (COMPLETE)

**Goal:** To independently verify the performance claims made in the ADMET-AI paper by running the provided tool on a public benchmark dataset.

**Outcome:** Phase 1 has been successfully completed. We have successfully reproduced the reported performance of the ADMET-AI model on the ClinTox benchmark dataset.

### Key Activities & Detailed Results:

1.  **Benchmark Identification:**
    *   **Decision:** We selected the **ClinTox** dataset from the Therapeutics Data Commons (TDC).
    *   **Reasoning:** This is a well-defined binary classification task mentioned in the paper's supplementary materials and is a good first test of the model's predictive power.

2.  **Data Download:**
    *   **Action:** A new directory `data/` was created to store all datasets.
    *   **Script:** A Python script, `src/download_data.py`, was written to download the data using the `PyTDC` library.
    *   **Output:** The script downloaded the ClinTox dataset and saved it to `data/data.csv`.

3.  **Prediction Execution:**
    *   **Initial Challenge:** The first attempt to run the `admet_predict` tool failed because the default SMILES column name ('smiles') was not present in the dataset.
    *   **Investigation:** We inspected the header of `data/data.csv` and found the SMILES data is in a column named **`Drug`**.
    *   **Command:** The prediction was successfully executed with the following command:
        ```bash
        admet_predict --data_path data/data.csv --smiles_column Drug --save_path data/predictions.csv
        ```
    *   **Output:** The predictions were saved to `data/predictions.csv`.

4.  **Results Analysis:**
    *   **Script:** A new Python script, `src/analyze_results.py`, was created to calculate the performance metric.
    *   **Methodology:** The script loads the true labels (the **`Y`** column) from `data/data.csv` and the predicted probabilities (the **`ClinTox`** column) from `data/predictions.csv`. It then uses the `scikit-learn` library to calculate the Area Under the Receiver Operating Characteristic Curve (AUROC).
    *   **Execution:** The script was run using `python3 src/analyze_results.py`.

5.  **Comparison & Conclusion:**
    *   **Our Calculated AUROC:** **0.9774**
    *   **Reported AUROC:** After locating the supplementary data (`btae416_supplementary_data/ADMET-AI Supplement.pdf`), we found the reported AUROC for the single-task ClinTox model in **Supplementary Figure 2, Panel C**. The value reported in the graph is approximately **0.98**.
    *   **Conclusion:** Our calculated AUROC score is in excellent agreement with the value reported in the paper. This successfully validates the paper's claims for this specific benchmark and concludes Phase 1.

---

## Phase 2: Novel Extension (Current Phase)

We are now beginning Phase 2, where we will build a novel extension to the ADMET-AI tool. The potential projects are:

1.  **Error Analysis & Visualization**
2.  **Comparative Drug Analysis**
3.  **New Dataset Testing**