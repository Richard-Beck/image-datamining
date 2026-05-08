#!/usr/bin/env bash
#SBATCH --job-name=k00_cpsam_dist
#SBATCH --array=1-320
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --qos=small
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
ANALYSIS_DIR="${REPO_ROOT}/analyses/K00_GemcitabineExposure_033023"
TASK_SCRIPT="${REPO_ROOT}/scripts/k00_cpsam_nearest_distances_stack.py"

MASK_DIR=${MASK_DIR:-"${ANALYSIS_DIR}/cpsam_full_stacks/masks"}
OUT_DIR=${OUT_DIR:-"${ANALYSIS_DIR}/cpsam_full_stacks/nearest_distances"}
SITE_MANIFEST=${SITE_MANIFEST:-"${OUT_DIR}/site_manifest.tsv"}
LOG_DIR=${LOG_DIR:-"${ANALYSIS_DIR}/logs/cpsam_nearest_distances"}
ARRAY_INDEX=${FILE_INDEX:-${SLURM_ARRAY_TASK_ID:-1}}

CONNECTIVITY=${CONNECTIVITY:-8}
START_FRAME=${START_FRAME:-0}
MAX_FRAMES=${MAX_FRAMES:-}
OVERWRITE=${OVERWRITE:-FALSE}

mkdir -p "$OUT_DIR" "$LOG_DIR"

if [ ! -s "$SITE_MANIFEST" ]; then
  MANIFEST_LOCK="${SITE_MANIFEST}.create.lock"
  if mkdir "$MANIFEST_LOCK" 2>/dev/null; then
    trap 'rm -rf "$MANIFEST_LOCK"' EXIT
    if [ ! -s "$SITE_MANIFEST" ]; then
      tmp_manifest="${SITE_MANIFEST}.tmp.$$"
      {
        printf "file_index\tsite_id\tmask_tiff\n"
        idx=0
        while IFS= read -r mask_tiff; do
          base=$(basename "$mask_tiff")
          site_id=${base%_cpsam_masks.tiff}
          site_id=${site_id%_cpsam_masks.tif}
          idx=$((idx + 1))
          printf "%s\t%s\t%s\n" "$idx" "$site_id" "$mask_tiff"
        done < <(find "$MASK_DIR" -maxdepth 1 -type f \( -name "*_cpsam_masks.tiff" -o -name "*_cpsam_masks.tif" \) | sort)
      } > "$tmp_manifest"
      mv "$tmp_manifest" "$SITE_MANIFEST"
    fi
    rm -rf "$MANIFEST_LOCK"
    trap - EXIT
  else
    while [ ! -s "$SITE_MANIFEST" ]; do
      sleep 2
    done
  fi
fi

MANIFEST_ROWS=$(awk 'END {print NR - 1}' "$SITE_MANIFEST")
SITE_ID=$(awk -F '\t' -v idx="$ARRAY_INDEX" 'NR == idx + 1 {print $2}' "$SITE_MANIFEST")
MASK_TIFF=$(awk -F '\t' -v idx="$ARRAY_INDEX" 'NR == idx + 1 {print $3}' "$SITE_MANIFEST")

TASK_LOG="${LOG_DIR}/site_${ARRAY_INDEX}_${SITE_ID:-unknown}.log"
exec > "$TASK_LOG" 2>&1

echo "[$(date --iso-8601=seconds)] Starting k00_cpsam_dist task"
echo "HOSTNAME=${HOSTNAME:-unknown}"
echo "PWD=$(pwd)"
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
echo "TASK_SCRIPT=${TASK_SCRIPT}"
echo "MASK_DIR=${MASK_DIR}"
echo "OUT_DIR=${OUT_DIR}"
echo "SITE_MANIFEST=${SITE_MANIFEST}"
echo "MANIFEST_ROWS=${MANIFEST_ROWS}"
echo "ARRAY_INDEX=${ARRAY_INDEX}"
echo "SITE_ID=${SITE_ID:-}"
echo "MASK_TIFF=${MASK_TIFF:-}"
echo "CONNECTIVITY=${CONNECTIVITY}"
echo "START_FRAME=${START_FRAME}"
echo "MAX_FRAMES=${MAX_FRAMES}"
echo "OVERWRITE=${OVERWRITE}"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"
echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID:-NA}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-NA}"
trap 'status=$?; echo "[$(date --iso-8601=seconds)] Task failed with exit status ${status} at line ${LINENO}"; exit ${status}' ERR

if [ "$ARRAY_INDEX" -gt "$MANIFEST_ROWS" ]; then
  echo "ARRAY_INDEX ${ARRAY_INDEX} is beyond ${MANIFEST_ROWS} manifest rows; exiting without work."
  exit 0
fi

if [ -z "$MASK_TIFF" ]; then
  echo "No mask TIFF found for ARRAY_INDEX ${ARRAY_INDEX}."
  exit 1
fi

CMD=(
  conda run -n cpose python "$TASK_SCRIPT"
  --mask-tiff "$MASK_TIFF"
  --out-dir "$OUT_DIR"
  --connectivity "$CONNECTIVITY"
  --start-frame "$START_FRAME"
)

if [ -n "$MAX_FRAMES" ]; then
  CMD+=(--max-frames "$MAX_FRAMES")
fi

if [ "$OVERWRITE" = "TRUE" ]; then
  CMD+=(--overwrite)
fi

echo "Running: ${CMD[*]}"
"${CMD[@]}"

echo "[$(date --iso-8601=seconds)] Finished k00_cpsam_dist task"
