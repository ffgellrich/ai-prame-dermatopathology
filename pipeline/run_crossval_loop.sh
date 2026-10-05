#!/usr/bin/env bash
# =============================================================================
# Project:     AI PRAME Study - Morphological Features & Malignancy
# Script:      run_crossval_loop.sh
# Author:      Frank Friedrich Gellrich
# Description: Executes 5-fold cross-validation across all diagnostic,
#              morphological, and immunohistochemical targets on the 
#              historical training cohort (80% split).
# =============================================================================

#SBATCH --job-name=stamp_cv80
#SBATCH --output=logs/slurm_cv80_%j.out
#SBATCH --error=logs/slurm_cv80_%j.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --time=120:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=20
#SBATCH --constraint=gpu24
set -euo pipefail

# Ensure central logs directory exists
mkdir -p logs

# Activate Python virtual environment
source .venv/bin/activate

# 1) Targets: 22 morphological criteria, diagnostic endpoints, and IHC markers
LABELS=(
  "KI_korr_Diagnose"
  "KI_maligne"
  "Junktionale_Beteiligung"
  "PRAME_Status"
  "Dermale_Beteiligung"
  "Symmetrie"
  "Aufsteigende_Melanozyten"
  "Pagetoide_Durchsetzung"
  "Konfluenz_junktionaler_Nester"
  "Ausgepraegte_junktionale_Atypien"
  "Atypische_junktionale_Mitosen"
  "Fehlende_Ausreifung_in_Nestern"
  "Ausgepraegte_dermale_Atypien"
  "Dermale_Mitosen"
  "Geordnete_Vereinzelung_Benignes_Muster"
  "Atypische_Vereinzelung_Dissoziation_Malignes_Muster"
  "Regression"
  "Dichtes_lymphozytaeres_Infiltrat"
  "Intratumorales_lymphozytaeres_Infiltrat_TILs"
  "Ulzeration"
  "Infiltration_der_Subkutis"
  "Mib1"
)

# 2) Paths relative to repository root
BASE_OUT="results/models/cross_val_HE_80"
CONFIG="config/config.yaml"

# 3) Cross-validation execution loop
for label in "${LABELS[@]}"; do
  outdir="${BASE_OUT}/${label}"
  
  # Skip feature if 5-fold cross-validation has already completed
  if [ -f "${outdir}/split-4/patient-preds.csv" ]; then
    echo "Skipping ${label}: output already exists at ${outdir}"
    continue
  fi

  mkdir -p "${outdir}"
  log_out="${outdir}/cross_val.${SLURM_JOB_ID:-local}_${label}.out"
  log_err="${outdir}/cross_val.${SLURM_JOB_ID:-local}_${label}.err"
  safe_label="${label// /_}"
  tmp_cfg="$(mktemp "/tmp/config.${safe_label}.XXXXXX")"

  # Dynamically inject label and output directory into temporary YAML configuration
  python - <<'PY' "${CONFIG}" "${tmp_cfg}" "${label}" "${outdir}"
import sys, yaml
src, dst, label, outdir = sys.argv[1:5]
with open(src, "r") as f:
    cfg = yaml.safe_load(f)
cfg.setdefault("crossval", {})
cfg["crossval"]["ground_truth_label"] = label
cfg["crossval"]["output_dir"] = outdir
with open(dst, "w") as f:
    yaml.safe_dump(cfg, f, sort_keys=False)
PY

  echo "Running 5-fold cross-validation for ${label} -> ${outdir}"
  stamp -c "${tmp_cfg}" crossval >"${log_out}" 2>"${log_err}"

  rm -f "${tmp_cfg}"
done

echo "Cross-validation loop completed successfully."