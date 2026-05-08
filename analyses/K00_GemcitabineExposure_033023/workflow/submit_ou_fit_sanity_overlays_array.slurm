#!/usr/bin/env bash
#SBATCH --job-name=k00_ou_gifs
#SBATCH --array=1-320
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --qos=small
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
ANALYSIS_DIR="${REPO_ROOT}/analyses/K00_GemcitabineExposure_033023"
TASK_SCRIPT="${ANALYSIS_DIR}/workflow/render_ou_fit_sanity_overlays.R"

SEED=${SEED:-17}
TRACKS_PER_GROUP=${TRACKS_PER_GROUP:-16}
CROP_BUFFER_PX=${CROP_BUFFER_PX:-40}
FPS=${FPS:-2}
ARRAY_INDEX=${FILE_INDEX:-${SLURM_ARRAY_TASK_ID:-1}}
QUIET=${QUIET:-FALSE}

INPUT_RDS=${INPUT_RDS:-"${ANALYSIS_DIR}/data/tracking_data_yellow_reconstructed_area2x_nonnegative_trackids_min3.rds"}
OUT_DIR=${OUT_DIR:-"${ANALYSIS_DIR}/overlay_checks/ou_fit_sanity/slurm_seed${SEED}_random_tracks_all_frames"}
RENDERED_CSV=${RENDERED_CSV:-"${OUT_DIR}/rendered_overlays.csv"}
LOG_DIR=${LOG_DIR:-"${ANALYSIS_DIR}/logs/ou_fit_sanity_overlays/seed${SEED}_random_tracks_all_frames"}

mkdir -p "$LOG_DIR"
TASK_LOG="${LOG_DIR}/track_${ARRAY_INDEX}.log"
exec > "$TASK_LOG" 2>&1

echo "[$(date --iso-8601=seconds)] Starting k00_ou_gifs task"
echo "HOSTNAME=${HOSTNAME:-unknown}"
echo "PWD=$(pwd)"
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
echo "INPUT_RDS=${INPUT_RDS}"
echo "OUT_DIR=${OUT_DIR}"
echo "RENDERED_CSV=${RENDERED_CSV}"
echo "SEED=${SEED}"
echo "TRACKS_PER_GROUP=${TRACKS_PER_GROUP}"
echo "CROP_BUFFER_PX=${CROP_BUFFER_PX}"
echo "FPS=${FPS}"
echo "ARRAY_INDEX=${ARRAY_INDEX}"
echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID:-NA}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-NA}"
trap 'status=$?; echo "[$(date --iso-8601=seconds)] Task failed with exit status ${status} at line ${LINENO}"; exit ${status}' ERR

"${REPO_ROOT}/scripts/agentRrunner.sh" "$TASK_SCRIPT" \
  --repo_root="${REPO_ROOT}" \
  --input_rds="${INPUT_RDS}" \
  --out_dir="${OUT_DIR}" \
  --rendered_csv="${RENDERED_CSV}" \
  --file_index="${ARRAY_INDEX}" \
  --tracks_per_group="${TRACKS_PER_GROUP}" \
  --crop_buffer_px="${CROP_BUFFER_PX}" \
  --seed="${SEED}" \
  --fps="${FPS}" \
  --quiet="${QUIET}"

echo "[$(date --iso-8601=seconds)] Finished k00_ou_gifs task"
