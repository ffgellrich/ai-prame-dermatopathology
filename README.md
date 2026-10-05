# AI PRAME Study: Artificial Intelligence for Dermatopathological Evaluation of Melanocytic Lesions

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.1.0-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-%E2%89%A5%203.9-blue.svg)](https://www.python.org/)

## Overview
This repository contains the complete code, configuration, and statistical analysis pipeline for the **AI PRAME Study**. The project evaluates the diagnostic performance of a foundation model (Virchow2) paired with a Vision Transformer (ViT) Multiple Instance Learning (MIL) architecture in predicting slide-level malignancy, specific morphological criteria, and immunohistochemical targets (PRAME, Ki-67) in melanocytic lesions.

The workflow is structured into a server-side AI pipeline (STAMP framework), a comprehensive R-based statistical analysis suite, and a Python-based visualization module for publication-ready figures.

## Repository Structure

```text
.
├── config/                  # AI hyperparameter settings and STAMP configurations
│   └── config.yaml          # Central pipeline config (Virchow2 extractor, ViT MIL parameters)
├── data/                    
│   ├── processed/           # Processed cross-validation and independent test set predictions
│   └── raw/                 # (Excluded from version control due to privacy constraints)
├── pipeline/                # Server-side bash scripts for SLURM execution
│   ├── feature_extraction.sh
│   ├── run_training_loop.sh
│   ├── run_crossval_loop.sh
│   └── run_deployment_loop.sh
├── scripts/                 # Core analysis and plotting scripts
│   ├── 01_Performance_Malignancy.R          # Slide-level ROC/PR analyses
│   ├── 03_SpiderPlots_AUC.R                 # Polar performance coordinate calculations
│   ├── 04_Combined_Model_Performance.R      # Morphology + Clinical metrics
│   ├── 05_Selected_Criteria_Model_Performance.R
│   ├── 06_Combined_ClinicPRAME_Model_Performance.R
│   ├── 07_Table1_Descriptive_Statistics.R   # Cohort demographics
│   ├── 08_Ground_Truth_Morphology_Benchmark.R
│   ├── 09_Frequency Comparison Table.R      # Human vs. AI frequency analysis
│   ├── 10_Supplementary_Four_Field_Analysis.R
│   └── generate_single_script.py            # Generates the 9-panel publication Master Grid (Figure 2)
└── results/                 # Generated outputs (fully reproducible from scripts)
    ├── csv/                 # Aggregated metrics and curve data
    ├── figures/             # High-resolution PNG and vector PDF plots
    └── tables/              # Formatted Word (.docx) tables for publication
```

## Prerequisites & Environment

### AI Pipeline (Server-side)
The model training and feature extraction rely on the **STAMP** (Spatial Transcriptomics & Artificial Intelligence in Pathology) framework.
* Python >= 3.9
* PyTorch >= 2.0 with CUDA support
* SLURM Workload Manager

### Statistical Analysis (Local)
**R Packages:**
* `tidyverse`, `tidymodels`, `pROC`, `PRROC`, `broom`
* `flextable`, `officer`, `gtsummary`
* `here` (for relative path resolution)

**Python Packages:**
* `pandas`, `numpy`, `matplotlib`, `seaborn`

## Execution Workflow

The complete analysis can be reproduced in three sequential stages:

### 1. AI Pipeline Execution (Server)
Run the SLURM bash scripts located in the `pipeline/` directory sequentially. These scripts utilize the parameters defined in `config/config.yaml` to extract 256 µm tiles using Virchow2, train the models (80% CV), and deploy inferences on the independent test set (20%). 
*(Note: Requires access to the raw WSI data, which is omitted from this public repository.)*

### 2. Statistical Analysis & Table Generation (R)
Execute the R scripts in the `scripts/` folder in numerical order (`01` through `10`). 
The scripts will read the raw prediction outputs from `data/processed/`, calculate all confidence intervals, ROC/PR curves, and logistic calibrations, and populate the `results/csv/` and `results/tables/` directories.

### 3. Master Figure Generation (Python)
Execute the final visualization script to generate the comprehensive 9-panel publication figure.
```bash
python scripts/generate_single_script.py
```
This script ingests the aggregated data from `results/csv/` and outputs high-resolution vector and raster graphics to `results/figures/Figure_2_Master_Grid.pdf`.

## Data Availability
Due to patient privacy regulations, the raw whole-slide images (WSI) and sensitive unanonymized clinical records are not publicly hosted. The `data/processed/` directory contains the fully anonymized, aggregated prediction scores necessary to reproduce all statistical tests, tables, and figures presented in the study.

## Author & Affiliation
**Frank Friedrich Gellrich, MD**  
Department of Dermatology, University Hospital Carl Gustav Carus Dresden  
Faculty of Medicine, TUD Dresden University of Technology  
Dresden, Germany

## License
This project is licensed under the MIT License - see the LICENSE file for details.