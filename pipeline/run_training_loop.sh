#!/usr/bin/env bash
# =============================================================================
# Project:     AI PRAME Study - Morphological Features & Malignancy
# Script:      run_training_loop.sh
# Author:      Frank Friedrich Gellrich
# Description: Executes sequential model training across all 22 diagnostic,
#              morphological, and immunohistochemical targets on the 
#              historical training cohort (80% split).
# =============================================================================

#SBATCH --job-name=stamp_training_loop
#SBATCH --output=logs/slurm_train_loop_%j.out
#SBATCH --error=logs/slurm_train_loop_%j.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --time=50:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=20
#SBATCH --constraint=gpu48
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
BASE_OUT="results/models/train_80_HE"
CONFIG="config/config.yaml"
FAILED_LABELS=()

# 3) Sequential training execution loop
for label in "${LABELS[@]}"; do
  outdir="${BASE_OUT}/${label}"
  
  # Skip feature if completed checkpoint already exists
  if [ -f "${outdir}/model.ckpt" ]; then
    echo "Skipping ${label}: completed model already exists at ${outdir}/model.ckpt"
    continue
  fi

  mkdir -p "${outdir}"
  log_out="${outdir}/training.${SLURM_JOB_ID:-local}_${label}.out"
  log_err="${outdir}/training.${SLURM_JOB_ID:-local}_${label}.err"
  safe_label="${label// /_}"
  tmp_cfg="$(mktemp "/tmp/config.${safe_label}.XXXXXX")"

  # Dynamically configure training target and output directory in temporary YAML
  python - <<'PY' "${CONFIG}" "${tmp_cfg}" "${label}" "${outdir}"
import sys, yaml
src, dst, label, outdir = sys.argv[1:5]
with open(src, "r") as f:
    cfg = yaml.safe_load(f)
cfg.setdefault("training", {})
cfg["training"]["ground_truth_label"] = label
cfg["training"]["output_dir"] = outdir
with open(dst, "w") as f:
    yaml.safe_dump(cfg, f, sort_keys=False)
PY

  echo "Running training for ${label} -> ${outdir}"
  echo "Per-label stdout: ${log_out}"
  echo "Per-label stderr: ${log_err}"

  # Run STAMP training with live console mirroring and per-feature file logging
  if stamp -c "${tmp_cfg}" train \
    > >(tee -a "${log_out}") \
    2> >(tee -a "${log_err}" >&2); then
    echo "Training completed for ${label}"
  else
    rc=$?
    echo "Training failed for ${label} with exit code ${rc}" >&2
    FAILED_LABELS+=("${label}")
  fi

  rm -f "${tmp_cfg}"
done

# Error handling and summary
if [ "${#FAILED_LABELS[@]}" -gt 0 ]; then
  echo "Training failed for labels: ${FAILED_LABELS[*]}" >&2
  exit 1
fi

echo "Training loop completed successfully."