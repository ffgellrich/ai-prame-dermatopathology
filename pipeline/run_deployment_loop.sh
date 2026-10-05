#!/usr/bin/env bash
# =============================================================================
# Project:     AI PRAME Study - Morphological Features & Malignancy
# Script:      run_deployment_loop.sh
# Author:      Frank Friedrich Gellrich
# Description: Executes model inference/deployment on the independent test
#              cohort (20% split) across all 22 diagnostic, morphological,
#              and immunohistochemical models trained on the 80% cohort.
# =============================================================================

#SBATCH --job-name=stamp_deploy20
#SBATCH --output=logs/slurm_deploy_%j.out
#SBATCH --error=logs/slurm_deploy_%j.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --time=10:00:00
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
BASE_OUT="results/models/deploy_20_HE"
CHECKPOINT="results/models/train_80_HE"
CONFIG="config/config.yaml"

# 3) Deployment execution loop
for label in "${LABELS[@]}"; do
  outdir="${BASE_OUT}/${label}"
  
  # Skip feature if deployment output folder already exists
  if [ -d "${outdir}" ]; then
    echo "Skipping ${label}: output dir already exists at ${outdir}"
    continue
  fi

  mkdir -p "${outdir}"
  checkpointdir="${CHECKPOINT}/${label}/model.ckpt"
  log_out="${outdir}/deploying.${SLURM_JOB_ID:-local}_${label}.out"
  log_err="${outdir}/deploying.${SLURM_JOB_ID:-local}_${label}.err"
  safe_label="${label// /_}"
  tmp_cfg="$(mktemp "/tmp/config.${safe_label}.XXXXXX")"

  # Dynamically configure deployment section in temporary YAML file
  python - <<'PY' "${CONFIG}" "${tmp_cfg}" "${label}" "${outdir}" "${checkpointdir}"
import sys, yaml
src, dst, label, outdir, checkpointdir = sys.argv[1:6]
with open(src, "r") as f:
    cfg = yaml.safe_load(f)
cfg.setdefault("deployment", {})
cfg["deployment"]["ground_truth_label"] = label
cfg["deployment"]["output_dir"] = outdir
cfg["deployment"]["checkpoint_paths"] = [checkpointdir]
with open(dst, "w") as f:
    yaml.safe_dump(cfg, f, sort_keys=False)
PY

  echo "Running deployment for ${label} -> ${outdir}"
  stamp -c "${tmp_cfg}" deploy >"${log_out}" 2>"${log_err}"

  rm -f "${tmp_cfg}"
done

echo "Deployment loop completed successfully."