#!/usr/bin/env bash
#SBATCH --job-name=k00_motility_censor
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --qos=small
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
ANALYSIS_DIR="${REPO_ROOT}/analyses/K00_GemcitabineExposure_033023"
LOG_DIR=${LOG_DIR:-"${ANALYSIS_DIR}/logs/motility_censoring_summaries"}
AT_RISK_RDS=${AT_RISK_RDS:-"${ANALYSIS_DIR}/data/yellow_censoring_at_risk.rds"}
MODELS_RDS=${MODELS_RDS:-"${ANALYSIS_DIR}/data/yellow_censoring_models.rds"}
OUT_DIR=${OUT_DIR:-"${ANALYSIS_DIR}/data/motility_censoring_summaries"}
MAX_LAG=${MAX_LAG:-5}
ACF_LAG=${ACF_LAG:-2}
FIXED_WINDOW=${FIXED_WINDOW:-5}
MAX_CENSOR_ROWS=${MAX_CENSOR_ROWS:-250000}
MAX_IPCW_WEIGHT=${MAX_IPCW_WEIGHT:-20}
SEED=${SEED:-20260509}

mkdir -p "$LOG_DIR"
TASK_LOG="${LOG_DIR}/motility_censoring_summaries_${SLURM_JOB_ID:-manual}.log"
exec > "$TASK_LOG" 2>&1

echo "[$(date --iso-8601=seconds)] Starting k00_motility_censor"
echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
echo "AT_RISK_RDS=${AT_RISK_RDS}"
echo "MODELS_RDS=${MODELS_RDS}"
echo "OUT_DIR=${OUT_DIR}"
trap 'status=$?; echo "[$(date --iso-8601=seconds)] Task failed with exit status ${status} at line ${LINENO}"; exit ${status}' ERR

"${REPO_ROOT}/scripts/agentRrunner.sh" "${ANALYSIS_DIR}/workflow/build_motility_censoring_summaries.R" \
  --analysis_dir="${ANALYSIS_DIR}" \
  --at_risk_rds="${AT_RISK_RDS}" \
  --censoring_model_rds="${MODELS_RDS}" \
  --out_dir="${OUT_DIR}" \
  --max_lag="${MAX_LAG}" \
  --acf_lag="${ACF_LAG}" \
  --fixed_window="${FIXED_WINDOW}" \
  --max_censor_rows="${MAX_CENSOR_ROWS}" \
  --max_ipcw_weight="${MAX_IPCW_WEIGHT}" \
  --seed="${SEED}"

echo "[$(date --iso-8601=seconds)] Finished k00_motility_censor"
