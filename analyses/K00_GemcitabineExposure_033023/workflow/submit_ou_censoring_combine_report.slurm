#!/usr/bin/env bash
#SBATCH --job-name=k00_ou_censor_report
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --qos=small
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
ANALYSIS_DIR="${REPO_ROOT}/analyses/K00_GemcitabineExposure_033023"
LOG_DIR=${LOG_DIR:-"${ANALYSIS_DIR}/logs/ou_censoring_combine_report"}
MODELS_RDS=${MODELS_RDS:-"${ANALYSIS_DIR}/data/yellow_censoring_models.rds"}
CONDITION_DIR=${CONDITION_DIR:-"${ANALYSIS_DIR}/data/ou_censoring_simulation_conditions"}
SIMULATION_RDS=${SIMULATION_RDS:-"${ANALYSIS_DIR}/data/ou_censoring_simulation.rds"}
MAX_LAG=${MAX_LAG:-5}
PARALLEL_CORES=${PARALLEL_CORES:-${SLURM_CPUS_PER_TASK:-32}}

mkdir -p "$LOG_DIR" "$(dirname "$SIMULATION_RDS")"
TASK_LOG="${LOG_DIR}/ou_censoring_combine_report_${SLURM_JOB_ID:-manual}.log"
exec > "$TASK_LOG" 2>&1

echo "[$(date --iso-8601=seconds)] Starting k00_ou_censor_report"
echo "HOSTNAME=${HOSTNAME:-unknown}"
echo "PWD=$(pwd)"
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
echo "MODELS_RDS=${MODELS_RDS}"
echo "CONDITION_DIR=${CONDITION_DIR}"
echo "SIMULATION_RDS=${SIMULATION_RDS}"
echo "MAX_LAG=${MAX_LAG}"
echo "PARALLEL_CORES=${PARALLEL_CORES}"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"
trap 'status=$?; echo "[$(date --iso-8601=seconds)] Task failed with exit status ${status} at line ${LINENO}"; exit ${status}' ERR

"${REPO_ROOT}/scripts/agentRrunner.sh" "${ANALYSIS_DIR}/workflow/combine_ou_censoring_conditions.R" \
  --analysis_dir="${ANALYSIS_DIR}" \
  --models_rds="${MODELS_RDS}" \
  --condition_dir="${CONDITION_DIR}" \
  --out_rds="${SIMULATION_RDS}" \
  --max_lag="${MAX_LAG}" \
  --parallel_cores="${PARALLEL_CORES}"

"${REPO_ROOT}/scripts/agentRrunner.sh" "${ANALYSIS_DIR}/workflow/render_yellow_censoring_simulation_report.R"

echo "[$(date --iso-8601=seconds)] Finished k00_ou_censor_report"
