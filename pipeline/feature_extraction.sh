#!/bin/bash
# =============================================================================
# Pipeline: Feature Extraction / Preprocessing (Virchow2)
# Framework: STAMP
# =============================================================================

#SBATCH --job-name=stamp_preprocess
#SBATCH --output=logs/slurm_preprocess_%j.out
#SBATCH --error=logs/slurm_preprocess_%j.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --time=150:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=40
#SBATCH --constraint=gpu24

# Ensure log directory exists
mkdir -p logs

# Activate Python virtual environment
# Adjust the path to your local virtual environment or conda installation
source .venv/bin/activate

# Execute STAMP feature extraction using config/config.yaml
# (Extracts tiles and generates Virchow2 embeddings)
stamp preprocess --config config/config.yaml