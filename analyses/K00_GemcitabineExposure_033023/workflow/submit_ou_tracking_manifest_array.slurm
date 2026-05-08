#!/usr/bin/env bash
#SBATCH --job-name=k00_ou_track
#SBATCH --array=1-500
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --qos=normal
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
ANALYSIS_DIR="${REPO_ROOT}/analyses/K00_GemcitabineExposure_033023"
MANIFEST_SCRIPT="${ANALYSIS_DIR}/workflow/make_ou_tracking_manifest.R"
TASK_SCRIPT="${ANALYSIS_DIR}/workflow/fit_ou_tracking_manifest_row.R"
LOG_DIR=${LOG_DIR:-"${ANALYSIS_DIR}/logs/ou_tracking_manifest_yellow_reconstructed_area2x_nonnegative_trackids_min3"}

INPUT_RDS=${INPUT_RDS:-"${ANALYSIS_DIR}/data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds"}
MANIFEST=${MANIFEST:-"${ANALYSIS_DIR}/data/ou_tracking_manifest_yellow_reconstructed_area2x_nonnegative_trackids_min3.csv"}
OUT_CSV=${OUT_CSV:-"${ANALYSIS_DIR}/data/ou_tracking_fits_yellow_reconstructed_area2x_nonnegative_trackids_min3.csv"}
FRAME_INTERVAL=${FRAME_INTERVAL:-2}
N_STARTS=${N_STARTS:-25}
ARRAY_INDEX=${FILE_INDEX:-${SLURM_ARRAY_TASK_ID:-1}}
QUIET=${QUIET:-FALSE}

if [ ! -s "$MANIFEST" ]; then
  MANIFEST_LOCK="${MANIFEST}.create.lock"
  if mkdir "$MANIFEST_LOCK" 2>/dev/null; then
    trap 'rm -rf "$MANIFEST_LOCK"' EXIT
    if [ ! -s "$MANIFEST" ]; then
      "${REPO_ROOT}/scripts/agentRrunner.sh" "$MANIFEST_SCRIPT" \
        --repo_root="${REPO_ROOT}" \
        --input_rds="${INPUT_RDS}" \
        --out_csv="${OUT_CSV}" \
        --out_manifest="${MANIFEST}" \
        --frame_interval="${FRAME_INTERVAL}" \
        --n_starts="${N_STARTS}" \
        --quiet="${QUIET}"
    fi
    rm -rf "$MANIFEST_LOCK"
    trap - EXIT
  else
    while [ ! -s "$MANIFEST" ]; do
      sleep 2
    done
  fi
fi

MANIFEST_ROWS=$(awk 'END {print NR - 1}' "$MANIFEST")
FILE_INDEX=$(( (ARRAY_INDEX - 1) / N_STARTS + 1 ))
START_ID=$(( (ARRAY_INDEX - 1) % N_STARTS + 1 ))

mkdir -p "$LOG_DIR"
TASK_LOG="${LOG_DIR}/condition_${FILE_INDEX}_start_${START_ID}.log"
exec > "$TASK_LOG" 2>&1

echo "[$(date --iso-8601=seconds)] Starting k00_ou_track task"
echo "HOSTNAME=${HOSTNAME:-unknown}"
echo "PWD=$(pwd)"
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
echo "INPUT_RDS=${INPUT_RDS}"
echo "MANIFEST=${MANIFEST}"
echo "MANIFEST_ROWS=${MANIFEST_ROWS}"
echo "OUT_CSV=${OUT_CSV}"
echo "FRAME_INTERVAL=${FRAME_INTERVAL}"
echo "N_STARTS=${N_STARTS}"
echo "ARRAY_INDEX=${ARRAY_INDEX}"
echo "FILE_INDEX=${FILE_INDEX}"
echo "START_ID=${START_ID}"
echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID:-NA}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-NA}"
trap 'status=$?; echo "[$(date --iso-8601=seconds)] Task failed with exit status ${status} at line ${LINENO}"; exit ${status}' ERR

if [ "$FILE_INDEX" -gt "$MANIFEST_ROWS" ]; then
  echo "FILE_INDEX ${FILE_INDEX} is beyond ${MANIFEST_ROWS} manifest rows; exiting without work."
  exit 0
fi

echo "Running fitter"
"${REPO_ROOT}/scripts/agentRrunner.sh" "$TASK_SCRIPT" \
  --repo_root="${REPO_ROOT}" \
  --manifest="${MANIFEST}" \
  --file_index="${FILE_INDEX}" \
  --n_starts="${N_STARTS}" \
  --start_id="${START_ID}" \
  --log_file="${TASK_LOG}" \
  --quiet="${QUIET}"

echo "[$(date --iso-8601=seconds)] Finished k00_ou_track task"
echo "Output CSV status:"
ls -lh "$OUT_CSV" 2>/dev/null || true
