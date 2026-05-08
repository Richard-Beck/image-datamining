#!/usr/bin/env bash
#SBATCH --job-name=k00_ou_grid
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --qos=small
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
ANALYSIS_DIR="${REPO_ROOT}/analyses/K00_GemcitabineExposure_033023"
TASK_SCRIPT="${ANALYSIS_DIR}/workflow/combine_ou_fit_sanity_overlay_grids.R"

SEED=${SEED:-17}
TRACKS_PER_PLOIDY=${TRACKS_PER_PLOIDY:-16}
GRID_COLS=${GRID_COLS:-4}
FPS=${FPS:-2}
ARRAY_INDEX=${FILE_INDEX:-${SLURM_ARRAY_TASK_ID:-1}}
QUIET=${QUIET:-FALSE}

OUT_DIR=${OUT_DIR:-"${ANALYSIS_DIR}/overlay_checks/ou_fit_sanity/slurm_seed${SEED}_random_tracks_all_frames"}
RENDERED_CSV=${RENDERED_CSV:-"${OUT_DIR}/rendered_overlays.csv"}
COMBINED_DIR=${COMBINED_DIR:-"${OUT_DIR}/combined"}
LOG_DIR=${LOG_DIR:-"${ANALYSIS_DIR}/logs/ou_fit_sanity_combine/seed${SEED}_random_tracks_all_frames"}

mkdir -p "$LOG_DIR"
TASK_LOG="${LOG_DIR}/dose_${ARRAY_INDEX}.log"
exec > "$TASK_LOG" 2>&1

echo "[$(date --iso-8601=seconds)] Starting k00_ou_grid task"
echo "HOSTNAME=${HOSTNAME:-unknown}"
echo "PWD=$(pwd)"
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
echo "OUT_DIR=${OUT_DIR}"
echo "RENDERED_CSV=${RENDERED_CSV}"
echo "COMBINED_DIR=${COMBINED_DIR}"
echo "TRACKS_PER_PLOIDY=${TRACKS_PER_PLOIDY}"
echo "GRID_COLS=${GRID_COLS}"
echo "FPS=${FPS}"
echo "ARRAY_INDEX=${ARRAY_INDEX}"
echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID:-NA}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-NA}"
trap 'status=$?; echo "[$(date --iso-8601=seconds)] Task failed with exit status ${status} at line ${LINENO}"; exit ${status}' ERR

"${REPO_ROOT}/scripts/agentRrunner.sh" "$TASK_SCRIPT" \
  --repo_root="${REPO_ROOT}" \
  --out_dir="${OUT_DIR}" \
  --rendered_csv="${RENDERED_CSV}" \
  --combined_dir="${COMBINED_DIR}" \
  --file_index="${ARRAY_INDEX}" \
  --tracks_per_ploidy="${TRACKS_PER_PLOIDY}" \
  --grid_cols="${GRID_COLS}" \
  --fps="${FPS}" \
  --quiet="${QUIET}"

echo "[$(date --iso-8601=seconds)] Finished k00_ou_grid task"
