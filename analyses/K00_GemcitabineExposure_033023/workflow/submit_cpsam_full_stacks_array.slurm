#!/usr/bin/env bash
#SBATCH --job-name=k00_cpsam
#SBATCH --array=1-320
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --qos=small
#SBATCH --gres=gpu:1
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
ANALYSIS_DIR="${REPO_ROOT}/analyses/K00_GemcitabineExposure_033023"
TASK_SCRIPT="${REPO_ROOT}/scripts/k00_run_cpsam_stack.py"

IMAGES_DIR=${IMAGES_DIR:-"/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/K00_GemcitabineExposure_033023/New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/Final_Tracking_analysis/Images_40Frames"}
OUT_DIR=${OUT_DIR:-"${ANALYSIS_DIR}/cpsam_full_stacks"}
SITE_MANIFEST=${SITE_MANIFEST:-"${OUT_DIR}/site_manifest.tsv"}
LOG_DIR=${LOG_DIR:-"${ANALYSIS_DIR}/logs/cpsam_full_stacks"}
ARRAY_INDEX=${FILE_INDEX:-${SLURM_ARRAY_TASK_ID:-1}}

CHANNEL_MODE=${CHANNEL_MODE:-rgb}
BATCH_SIZE=${BATCH_SIZE:-4}
MIN_SIZE=${MIN_SIZE:-15}
FLOW_THRESHOLD=${FLOW_THRESHOLD:-0.4}
CELLPROB_THRESHOLD=${CELLPROB_THRESHOLD:-0.0}
DIAMETER=${DIAMETER:-}
OVERWRITE=${OVERWRITE:-FALSE}

mkdir -p "$OUT_DIR" "$LOG_DIR"

if [ ! -s "$SITE_MANIFEST" ]; then
  MANIFEST_LOCK="${SITE_MANIFEST}.create.lock"
  if mkdir "$MANIFEST_LOCK" 2>/dev/null; then
    trap 'rm -rf "$MANIFEST_LOCK"' EXIT
    if [ ! -s "$SITE_MANIFEST" ]; then
      {
        printf "file_index\tinput_tiff\n"
        find "$IMAGES_DIR" -maxdepth 1 -type f \( -name "*.tiff" -o -name "*.tif" \) \
          | sort \
          | awk '{print NR "\t" $0}'
      } > "$SITE_MANIFEST"
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
INPUT_TIFF=$(awk -F '\t' -v idx="$ARRAY_INDEX" 'NR == idx + 1 {print $2}' "$SITE_MANIFEST")

TASK_LOG="${LOG_DIR}/site_${ARRAY_INDEX}.log"
exec > "$TASK_LOG" 2>&1

echo "[$(date --iso-8601=seconds)] Starting k00_cpsam task"
echo "HOSTNAME=${HOSTNAME:-unknown}"
echo "PWD=$(pwd)"
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
echo "TASK_SCRIPT=${TASK_SCRIPT}"
echo "IMAGES_DIR=${IMAGES_DIR}"
echo "OUT_DIR=${OUT_DIR}"
echo "SITE_MANIFEST=${SITE_MANIFEST}"
echo "MANIFEST_ROWS=${MANIFEST_ROWS}"
echo "ARRAY_INDEX=${ARRAY_INDEX}"
echo "INPUT_TIFF=${INPUT_TIFF:-}"
echo "CHANNEL_MODE=${CHANNEL_MODE}"
echo "BATCH_SIZE=${BATCH_SIZE}"
echo "MIN_SIZE=${MIN_SIZE}"
echo "FLOW_THRESHOLD=${FLOW_THRESHOLD}"
echo "CELLPROB_THRESHOLD=${CELLPROB_THRESHOLD}"
echo "DIAMETER=${DIAMETER}"
echo "OVERWRITE=${OVERWRITE}"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"
echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID:-NA}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-NA}"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-NA}"
trap 'status=$?; echo "[$(date --iso-8601=seconds)] Task failed with exit status ${status} at line ${LINENO}"; exit ${status}' ERR

if [ "$ARRAY_INDEX" -gt "$MANIFEST_ROWS" ]; then
  echo "ARRAY_INDEX ${ARRAY_INDEX} is beyond ${MANIFEST_ROWS} manifest rows; exiting without work."
  exit 0
fi

if [ -z "$INPUT_TIFF" ]; then
  echo "No input TIFF found for ARRAY_INDEX ${ARRAY_INDEX}."
  exit 1
fi

CMD=(
  conda run -n cpose python "$TASK_SCRIPT"
  --input-tiff "$INPUT_TIFF"
  --out-dir "$OUT_DIR"
  --channel-mode "$CHANNEL_MODE"
  --batch-size "$BATCH_SIZE"
  --min-size "$MIN_SIZE"
  --flow-threshold "$FLOW_THRESHOLD"
  --cellprob-threshold "$CELLPROB_THRESHOLD"
)

if [ -n "$DIAMETER" ]; then
  CMD+=(--diameter "$DIAMETER")
fi

if [ "$OVERWRITE" = "TRUE" ]; then
  CMD+=(--overwrite)
fi

echo "Running: ${CMD[*]}"
"${CMD[@]}"

echo "[$(date --iso-8601=seconds)] Finished k00_cpsam task"
