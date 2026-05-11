#!/usr/bin/env bash
#SBATCH --job-name=k00_yellow_censor_flow
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --qos=small
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
ANALYSIS_DIR="${REPO_ROOT}/analyses/K00_GemcitabineExposure_033023"
LOG_DIR=${LOG_DIR:-"${ANALYSIS_DIR}/logs/yellow_censoring_motility_workflow"}

RUN_SITE_ARRAY=${RUN_SITE_ARRAY:-TRUE}
RUN_CENSORING_MODEL=${RUN_CENSORING_MODEL:-TRUE}
RUN_OU_SIMULATION=${RUN_OU_SIMULATION:-TRUE}
RUN_MOTILITY=${RUN_MOTILITY:-TRUE}
RENDER_REPORT=${RENDER_REPORT:-TRUE}

SITE_ROOT=${SITE_ROOT:-"${ANALYSIS_DIR}/data/yellow_censoring_sites"}
SITE_DIR=${SITE_DIR:-"${SITE_ROOT}/sites"}
AT_RISK_RDS=${AT_RISK_RDS:-"${ANALYSIS_DIR}/data/yellow_censoring_at_risk.rds"}
STAGE_SUMMARY_CSV=${STAGE_SUMMARY_CSV:-"${ANALYSIS_DIR}/data/yellow_censoring_stage_summary.csv"}
REASON_SUMMARY_CSV=${REASON_SUMMARY_CSV:-"${ANALYSIS_DIR}/data/yellow_censoring_reason_summary.csv"}
MODELS_RDS=${MODELS_RDS:-"${ANALYSIS_DIR}/data/yellow_censoring_models.rds"}
OU_CONDITION_DIR=${OU_CONDITION_DIR:-"${ANALYSIS_DIR}/data/ou_censoring_simulation_conditions"}
OU_SIMULATION_RDS=${OU_SIMULATION_RDS:-"${ANALYSIS_DIR}/data/ou_censoring_simulation.rds"}
MOTILITY_DIR=${MOTILITY_DIR:-"${ANALYSIS_DIR}/data/motility_censoring_summaries"}
MAX_LAG=${MAX_LAG:-5}
MAX_CENSOR_ROWS=${MAX_CENSOR_ROWS:-250000}
MODEL_MAX_ROWS=${MODEL_MAX_ROWS:-0}
SEED=${SEED:-20260509}

mkdir -p "$LOG_DIR"
TASK_LOG="${LOG_DIR}/yellow_censoring_motility_workflow_${SLURM_JOB_ID:-manual}.log"
exec > "$TASK_LOG" 2>&1

as_true() {
  case "${1:-}" in
    1|TRUE|true|True|yes|YES|y|Y) return 0 ;;
    *) return 1 ;;
  esac
}

require_file() {
  if [ ! -s "$1" ]; then
    echo "Required file is missing or empty: $1" >&2
    exit 1
  fi
}

echo "[$(date --iso-8601=seconds)] Starting k00_yellow_censor_flow"
echo "HOSTNAME=${HOSTNAME:-unknown}"
echo "PWD=$(pwd)"
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
echo "RUN_SITE_ARRAY=${RUN_SITE_ARRAY}"
echo "RUN_CENSORING_MODEL=${RUN_CENSORING_MODEL}"
echo "RUN_OU_SIMULATION=${RUN_OU_SIMULATION}"
echo "RUN_MOTILITY=${RUN_MOTILITY}"
echo "RENDER_REPORT=${RENDER_REPORT}"
echo "SITE_DIR=${SITE_DIR}"
echo "AT_RISK_RDS=${AT_RISK_RDS}"
echo "MODELS_RDS=${MODELS_RDS}"
echo "OU_CONDITION_DIR=${OU_CONDITION_DIR}"
echo "OU_SIMULATION_RDS=${OU_SIMULATION_RDS}"
echo "MOTILITY_DIR=${MOTILITY_DIR}"
echo "MAX_LAG=${MAX_LAG}"
echo "MAX_CENSOR_ROWS=${MAX_CENSOR_ROWS}"
echo "MODEL_MAX_ROWS=${MODEL_MAX_ROWS}"
echo "SEED=${SEED}"
trap 'status=$?; echo "[$(date --iso-8601=seconds)] Workflow failed with exit status ${status} at line ${LINENO}"; exit ${status}' ERR

if as_true "$RUN_SITE_ARRAY"; then
  echo "[$(date --iso-8601=seconds)] Submitting yellow censoring site array"
  sbatch --wait --export=ALL,REPO_ROOT="${REPO_ROOT}",STAGED_DIR="${SITE_ROOT}" \
    "${ANALYSIS_DIR}/workflow/submit_yellow_censoring_sites_array.slurm"
else
  echo "[$(date --iso-8601=seconds)] Skipping site array; using existing site artifacts"
fi

if as_true "$RUN_CENSORING_MODEL"; then
  "${REPO_ROOT}/scripts/agentRrunner.sh" "${ANALYSIS_DIR}/workflow/combine_yellow_censoring_sites.R" \
    --analysis_dir="${ANALYSIS_DIR}" \
    --site_dir="${SITE_DIR}" \
    --out_rds="${AT_RISK_RDS}" \
    --out_summary_csv="${STAGE_SUMMARY_CSV}" \
    --out_reason_csv="${REASON_SUMMARY_CSV}"

  "${REPO_ROOT}/scripts/agentRrunner.sh" "${ANALYSIS_DIR}/workflow/fit_yellow_censoring_models.R" \
    --analysis_dir="${ANALYSIS_DIR}" \
    --at_risk_rds="${AT_RISK_RDS}" \
    --out_rds="${MODELS_RDS}" \
    --max_rows="${MODEL_MAX_ROWS}" \
    --seed="${SEED}"
else
  echo "[$(date --iso-8601=seconds)] Skipping censoring combine/model; using existing artifacts"
  require_file "$AT_RISK_RDS"
  require_file "$MODELS_RDS"
fi

require_file "$MODELS_RDS"

if as_true "$RUN_OU_SIMULATION"; then
  echo "[$(date --iso-8601=seconds)] Submitting OU censoring condition array"
  sbatch --wait --export=ALL,REPO_ROOT="${REPO_ROOT}",CONDITION_DIR="${OU_CONDITION_DIR}",MODELS_RDS="${MODELS_RDS}",MAX_LAG="${MAX_LAG}",SEED="${SEED}" \
    "${ANALYSIS_DIR}/workflow/submit_ou_censoring_conditions_array.slurm"

  "${REPO_ROOT}/scripts/agentRrunner.sh" "${ANALYSIS_DIR}/workflow/combine_ou_censoring_conditions.R" \
    --analysis_dir="${ANALYSIS_DIR}" \
    --models_rds="${MODELS_RDS}" \
    --condition_dir="${OU_CONDITION_DIR}" \
    --out_rds="${OU_SIMULATION_RDS}" \
    --max_lag="${MAX_LAG}" \
    --parallel_cores="${SLURM_CPUS_PER_TASK:-1}"
else
  echo "[$(date --iso-8601=seconds)] Skipping OU simulation; using existing combined artifact"
fi

if as_true "$RUN_MOTILITY"; then
  "${REPO_ROOT}/scripts/agentRrunner.sh" "${ANALYSIS_DIR}/workflow/build_motility_censoring_summaries.R" \
    --analysis_dir="${ANALYSIS_DIR}" \
    --at_risk_rds="${AT_RISK_RDS}" \
    --censoring_model_rds="${MODELS_RDS}" \
    --out_dir="${MOTILITY_DIR}" \
    --max_lag="${MAX_LAG}" \
    --max_censor_rows="${MAX_CENSOR_ROWS}" \
    --seed="${SEED}"
fi

if as_true "$RENDER_REPORT"; then
  require_file "$OU_SIMULATION_RDS"
  "${REPO_ROOT}/scripts/agentRrunner.sh" "${ANALYSIS_DIR}/workflow/render_yellow_censoring_simulation_report.R"
fi

echo "[$(date --iso-8601=seconds)] Finished k00_yellow_censor_flow"
