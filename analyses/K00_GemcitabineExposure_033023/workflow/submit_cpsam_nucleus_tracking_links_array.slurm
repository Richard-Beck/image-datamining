#!/usr/bin/env bash
#SBATCH --job-name=k00_cpsam_link
#SBATCH --array=1-320
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --qos=small
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
ANALYSIS_DIR="${REPO_ROOT}/analyses/K00_GemcitabineExposure_033023"
TASK_SCRIPT="${ANALYSIS_DIR}/workflow/prototype_link_cpsam_to_nucleus_tracks.R"

RAW_TIFF_DIR=${RAW_TIFF_DIR:-"/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/K00_GemcitabineExposure_033023/New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/Final_Tracking_analysis/Images_40Frames"}
TRACKING_DIR=${TRACKING_DIR:-"/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/BreastCancerCellLines/SUM-159/K00_GemcitabineExposure_033023/New_20240125_SUM159_2N_4N_Gemcitabine_Incucyte_2hr(Analysis_QI_Core)/Final_Tracking_analysis/Tracking_CSVs"}
MASK_DIR=${MASK_DIR:-"${ANALYSIS_DIR}/cpsam_full_stacks/masks"}
OUT_DIR=${OUT_DIR:-"${ANALYSIS_DIR}/cpsam_nucleus_tracking_links"}
SITE_MANIFEST=${SITE_MANIFEST:-"${OUT_DIR}/site_manifest.tsv"}
LOG_DIR=${LOG_DIR:-"${ANALYSIS_DIR}/logs/cpsam_nucleus_tracking_links"}
ARRAY_INDEX=${FILE_INDEX:-${SLURM_ARRAY_TASK_ID:-1}}

FRAMES=${FRAMES:-all}
TOUCH_CONNECTIVITY=${TOUCH_CONNECTIVITY:-8}
OVERWRITE=${OVERWRITE:-FALSE}
TEXT_CEX=${TEXT_CEX:-0.95}

mkdir -p "$OUT_DIR" "$LOG_DIR"

if [ ! -s "$SITE_MANIFEST" ]; then
  MANIFEST_LOCK="${SITE_MANIFEST}.create.lock"
  if mkdir "$MANIFEST_LOCK" 2>/dev/null; then
    trap 'rm -rf "$MANIFEST_LOCK"' EXIT
    if [ ! -s "$SITE_MANIFEST" ]; then
      tmp_manifest="${SITE_MANIFEST}.tmp.$$"
      {
        printf "file_index\tsite_id\traw_tiff\tmask_tiff\ttracks_csv\n"
        idx=0
        while IFS= read -r mask_tiff; do
          base=$(basename "$mask_tiff")
          site_id=${base%_cpsam_masks.tiff}
          site_id=${site_id%_cpsam_masks.tif}
          raw_tiff="${RAW_TIFF_DIR}/${site_id}.tiff"
          if [ ! -f "$raw_tiff" ]; then
            raw_tiff="${RAW_TIFF_DIR}/${site_id}.tif"
          fi
          tracks_csv="${TRACKING_DIR}/${site_id}_CSV-Table.tiff.csv"
          if [ -f "$raw_tiff" ] && [ -f "$tracks_csv" ]; then
            idx=$((idx + 1))
            printf "%s\t%s\t%s\t%s\t%s\n" "$idx" "$site_id" "$raw_tiff" "$mask_tiff" "$tracks_csv"
          fi
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

TASK_LOG="${LOG_DIR}/site_${ARRAY_INDEX}_${SITE_ID:-unknown}.log"
exec > "$TASK_LOG" 2>&1

echo "[$(date --iso-8601=seconds)] Starting k00_cpsam_link task"
echo "HOSTNAME=${HOSTNAME:-unknown}"
echo "PWD=$(pwd)"
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
echo "TASK_SCRIPT=${TASK_SCRIPT}"
echo "RAW_TIFF_DIR=${RAW_TIFF_DIR}"
echo "TRACKING_DIR=${TRACKING_DIR}"
echo "MASK_DIR=${MASK_DIR}"
echo "OUT_DIR=${OUT_DIR}"
echo "SITE_MANIFEST=${SITE_MANIFEST}"
echo "MANIFEST_ROWS=${MANIFEST_ROWS}"
echo "ARRAY_INDEX=${ARRAY_INDEX}"
echo "SITE_ID=${SITE_ID:-}"
echo "FRAMES=${FRAMES}"
echo "TOUCH_CONNECTIVITY=${TOUCH_CONNECTIVITY}"
echo "OVERWRITE=${OVERWRITE}"
echo "TEXT_CEX=${TEXT_CEX}"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"
echo "SLURM_ARRAY_JOB_ID=${SLURM_ARRAY_JOB_ID:-NA}"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-NA}"
trap 'status=$?; echo "[$(date --iso-8601=seconds)] Task failed with exit status ${status} at line ${LINENO}"; exit ${status}' ERR

if [ "$ARRAY_INDEX" -gt "$MANIFEST_ROWS" ]; then
  echo "ARRAY_INDEX ${ARRAY_INDEX} is beyond ${MANIFEST_ROWS} manifest rows; exiting without work."
  exit 0
fi

if [ -z "$SITE_ID" ]; then
  echo "No site found for ARRAY_INDEX ${ARRAY_INDEX}."
  exit 1
fi

CMD=(
  scripts/agentRrunner.sh "$TASK_SCRIPT"
  --repo_root="${REPO_ROOT}"
  --manifest="${SITE_MANIFEST}"
  --file_index="${ARRAY_INDEX}"
  --out_dir="${OUT_DIR}"
  --frames="${FRAMES}"
  --touch_connectivity="${TOUCH_CONNECTIVITY}"
  --overwrite="${OVERWRITE}"
  --text_cex="${TEXT_CEX}"
)

echo "Running: ${CMD[*]}"
"${CMD[@]}"

echo "[$(date --iso-8601=seconds)] Finished k00_cpsam_link task"
