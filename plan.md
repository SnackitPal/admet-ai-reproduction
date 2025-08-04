# Strategic Plan: ADMET-AI Reproduction & Extension

This document outlines the phased approach for this project. The goal is to first reproduce the scientific claims of the ADMET-AI paper and then to build a novel extension.

---

### Phase 0: Environment Setup & Configuration (COMPLETE)

**Goal:** Establish a stable, correct, and reproducible development environment.

**Key Activities:**
*   Set up WSL/Ubuntu environment with native Node.js and Gemini CLI.
*   Established a Conda environment (`admet_project`) using a stable Python version (3.11).
*   Successfully resolved all dependency and compiler issues (`build-essential`, PyTorch version conflicts).
*   Installed the `admet-ai` package and verified its functionality with a "smoke test" (`admet_predict --help`).
*   Initialized a Git repository for version control.

---

### Phase 1: Baseline Scientific Reproduction (CURRENT PHASE)

**Goal:** To independently verify the performance claims made in the ADMET-AI paper by running the provided tool on a public benchmark dataset.

**Key Activities:**
1.  **Identify Benchmark:** Select a suitable benchmark dataset from the Therapeutics Data Commons (TDC) that was used in the original paper (e.g., ClinTox, hERG, etc.).
2.  **Download Data:** Fetch the public benchmark dataset onto the local machine.
3.  **Execute Prediction:** Run the `admet_predict` command-line tool on the entire dataset to generate predictions.
4.  **Analyze Results:** Write a new Python script (`src/analyze_results.py`) to calculate the relevant performance metric (e.g., AUC, Accuracy) from the prediction output.
5.  **Compare & Conclude:** Compare the calculated metric against the value reported in the paper's tables. A successful reproduction is one where the results are closely aligned.

---

### Phase 2: Novel Extension (FUTURE PHASE)

**Goal:** To build upon the reproduced work by adding a unique contribution. This demonstrates independent research capability.

**Potential Extension Projects (Choose ONE):**
1.  **Error Analysis & Visualization:** Identify the top 5 molecules that the model gets wrong (false positives/negatives). Write a script to generate and save images of their chemical structures using RDKit for analysis.
2.  **Comparative Drug Analysis:** Add a feature to compare a new molecule's predicted ADMET profile against a well-known drug (e.g., Metformin, Aspirin).
3.  **New Dataset Testing:** Test the pre-trained ADMET-AI model on a completely different, new public dataset of molecules it has never seen before to test its generalization capability.

---

### Phase 3: Documentation & Presentation (FINAL PHASE)

**Goal:** To package the project for professional presentation (e.g., on a resume, in an interview).

**Key Activities:**
1.  **Finalize `README.md`:** Create a high-quality, detailed README for the project's GitHub repository, explaining the project, the reproduction results, the extension, and how to run the code.
2.  **Write Project Article:** Author a blog post or a short article summarizing the project's findings and learning journey.

