#!/usr/bin/env bash
#SBATCH --job-name=k00_yellow_fit
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --qos=small
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
ANALYSIS_DIR="${REPO_ROOT}/analyses/K00_GemcitabineExposure_033023"
LOG_DIR=${LOG_DIR:-"${ANALYSIS_DIR}/logs/yellow_censoring_model_fit"}
SITE_ROOT=${SITE_ROOT:-"${ANALYSIS_DIR}/data/yellow_censoring_sites"}
SITE_DIR=${SITE_DIR:-"${SITE_ROOT}/sites"}
AT_RISK_RDS=${AT_RISK_RDS:-"${ANALYSIS_DIR}/data/yellow_censoring_at_risk.rds"}
STAGE_SUMMARY_CSV=${STAGE_SUMMARY_CSV:-"${ANALYSIS_DIR}/data/yellow_censoring_stage_summary.csv"}
REASON_SUMMARY_CSV=${REASON_SUMMARY_CSV:-"${ANALYSIS_DIR}/data/yellow_censoring_reason_summary.csv"}
MODEL_OUT=${MODEL_OUT:-"${ANALYSIS_DIR}/data/yellow_censoring_models.rds"}
MAX_ROWS=${MAX_ROWS:-0}
SEED=${SEED:-20260509}

mkdir -p "$LOG_DIR" "$(dirname "$MODEL_OUT")"
TASK_LOG="${LOG_DIR}/yellow_censoring_model_fit_${SLURM_JOB_ID:-manual}.log"
exec > "$TASK_LOG" 2>&1

echo "[$(date --iso-8601=seconds)] Starting k00_yellow_fit"
echo "HOSTNAME=${HOSTNAME:-unknown}"
echo "PWD=$(pwd)"
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
echo "SITE_DIR=${SITE_DIR}"
echo "AT_RISK_RDS=${AT_RISK_RDS}"
echo "MODEL_OUT=${MODEL_OUT}"
echo "MAX_ROWS=${MAX_ROWS}"
echo "SEED=${SEED}"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"
trap 'status=$?; echo "[$(date --iso-8601=seconds)] Task failed with exit status ${status} at line ${LINENO}"; exit ${status}' ERR

"${REPO_ROOT}/scripts/agentRrunner.sh" "${ANALYSIS_DIR}/workflow/combine_yellow_censoring_sites.R" \
  --analysis_dir="${ANALYSIS_DIR}" \
  --site_dir="${SITE_DIR}" \
  --out_rds="${AT_RISK_RDS}" \
  --out_summary_csv="${STAGE_SUMMARY_CSV}" \
  --out_reason_csv="${REASON_SUMMARY_CSV}"

"${REPO_ROOT}/scripts/agentRrunner.sh" "${ANALYSIS_DIR}/workflow/fit_yellow_censoring_models.R" \
  --analysis_dir="${ANALYSIS_DIR}" \
  --at_risk_rds="${AT_RISK_RDS}" \
  --out_rds="${MODEL_OUT}" \
  --max_rows="${MAX_ROWS}" \
  --seed="${SEED}"

echo "[$(date --iso-8601=seconds)] Finished k00_yellow_fit"
