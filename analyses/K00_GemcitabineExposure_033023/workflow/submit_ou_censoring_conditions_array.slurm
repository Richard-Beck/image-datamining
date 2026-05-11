#!/usr/bin/env bash
#SBATCH --job-name=k00_ou_censor
#SBATCH --array=1-20
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --qos=small
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
ANALYSIS_DIR="${REPO_ROOT}/analyses/K00_GemcitabineExposure_033023"
TASK_SCRIPT="${ANALYSIS_DIR}/workflow/simulate_ou_censoring_condition.R"

CONDITION_DIR=${CONDITION_DIR:-"${ANALYSIS_DIR}/data/ou_censoring_simulation_conditions"}
CONDITION_MANIFEST=${CONDITION_MANIFEST:-"${CONDITION_DIR}/condition_manifest.tsv"}
MODELS_RDS=${MODELS_RDS:-"${ANALYSIS_DIR}/data/yellow_censoring_models.rds"}
LOG_DIR=${LOG_DIR:-"${ANALYSIS_DIR}/logs/ou_censoring_conditions"}
ARRAY_INDEX=${FILE_INDEX:-${SLURM_ARRAY_TASK_ID:-1}}

BATCH_TRACKS=${BATCH_TRACKS:-1000}
PARALLEL_BATCHES=${PARALLEL_BATCHES:-${SLURM_CPUS_PER_TASK:-8}}
UNCENSORED_SAMPLE_TRACKS=${UNCENSORED_SAMPLE_TRACKS:-1000}
INCLUDE_CONTEXT_MODEL=${INCLUDE_CONTEXT_MODEL:-FALSE}
MIN_SEGMENT_FRAMES=${MIN_SEGMENT_FRAMES:-3}
MAX_FRAME=${MAX_FRAME:-39}
MAX_LAG=${MAX_LAG:-10}
SEED=${SEED:-20260509}
OVERWRITE=${OVERWRITE:-FALSE}

mkdir -p "$CONDITION_DIR" "$LOG_DIR"

if [ ! -s "$CONDITION_MANIFEST" ]; then
  tmp_manifest="${CONDITION_MANIFEST}.tmp.$$"
  {
    printf "file_index\tcondition\n"
    idx=0
    for ploidy in 2N 4N; do
      for dose in 0 3.125 6.25 12.5 25 50 100 200 400 800; do
        idx=$((idx + 1))
        printf "%s\t%s:%s\n" "$idx" "$ploidy" "$dose"
      done
    done
  } > "$tmp_manifest"
  mv "$tmp_manifest" "$CONDITION_MANIFEST"
fi

MANIFEST_ROWS=$(awk 'END {print NR - 1}' "$CONDITION_MANIFEST")
CONDITION=$(awk -F '\t' -v idx="$ARRAY_INDEX" 'NR == idx + 1 {print $2}' "$CONDITION_MANIFEST")

TASK_LOG="${LOG_DIR}/condition_${ARRAY_INDEX}_${CONDITION//:/_}.log"
exec > "$TASK_LOG" 2>&1

echo "[$(date --iso-8601=seconds)] Starting k00_ou_censor task"
echo "HOSTNAME=${HOSTNAME:-unknown}"
echo "PWD=$(pwd)"
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
echo "TASK_SCRIPT=${TASK_SCRIPT}"
echo "MODELS_RDS=${MODELS_RDS}"
echo "CONDITION_DIR=${CONDITION_DIR}"
echo "CONDITION_MANIFEST=${CONDITION_MANIFEST}"
echo "MANIFEST_ROWS=${MANIFEST_ROWS}"
echo "ARRAY_INDEX=${ARRAY_INDEX}"
echo "CONDITION=${CONDITION:-}"
echo "BATCH_TRACKS=${BATCH_TRACKS}"
echo "PARALLEL_BATCHES=${PARALLEL_BATCHES}"
echo "INCLUDE_CONTEXT_MODEL=${INCLUDE_CONTEXT_MODEL}"
echo "MIN_SEGMENT_FRAMES=${MIN_SEGMENT_FRAMES}"
echo "SEED=${SEED}"
echo "OVERWRITE=${OVERWRITE}"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"
echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID:-NA}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-NA}"
trap 'status=$?; echo "[$(date --iso-8601=seconds)] Task failed with exit status ${status} at line ${LINENO}"; exit ${status}' ERR

if [ "$ARRAY_INDEX" -gt "$MANIFEST_ROWS" ]; then
  echo "ARRAY_INDEX ${ARRAY_INDEX} is beyond ${MANIFEST_ROWS} manifest rows; exiting without work."
  exit 0
fi

"${REPO_ROOT}/scripts/agentRrunner.sh" "$TASK_SCRIPT" \
  --analysis_dir="${ANALYSIS_DIR}" \
  --models_rds="${MODELS_RDS}" \
  --manifest="${CONDITION_MANIFEST}" \
  --file_index="${ARRAY_INDEX}" \
  --out_dir="${CONDITION_DIR}" \
  --batch_tracks="${BATCH_TRACKS}" \
  --parallel_batches="${PARALLEL_BATCHES}" \
  --uncensored_sample_tracks="${UNCENSORED_SAMPLE_TRACKS}" \
  --include_context_model="${INCLUDE_CONTEXT_MODEL}" \
  --min_segment_frames="${MIN_SEGMENT_FRAMES}" \
  --max_frame="${MAX_FRAME}" \
  --max_lag="${MAX_LAG}" \
  --seed="${SEED}" \
  --overwrite="${OVERWRITE}"

echo "[$(date --iso-8601=seconds)] Finished k00_ou_censor task"
