#!/usr/bin/env bash
#SBATCH --job-name=k00_yellow_censor
#SBATCH --array=1-320
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --qos=small
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

set -euo pipefail

REPO_ROOT=${REPO_ROOT:-$(pwd)}
ANALYSIS_DIR="${REPO_ROOT}/analyses/K00_GemcitabineExposure_033023"
TASK_SCRIPT="${ANALYSIS_DIR}/workflow/build_yellow_censoring_site.R"

TRACKPOINT_LINKS_DIR=${TRACKPOINT_LINKS_DIR:-"${ANALYSIS_DIR}/cpsam_nucleus_tracking_links/trackpoint_links"}
NEAREST_DISTANCES_DIR=${NEAREST_DISTANCES_DIR:-"${ANALYSIS_DIR}/cpsam_full_stacks/nearest_distances/object_tables"}
STAGED_DIR=${STAGED_DIR:-"${ANALYSIS_DIR}/data/yellow_censoring_sites"}
SITE_DIR=${SITE_DIR:-"${STAGED_DIR}/sites"}
SITE_MANIFEST=${SITE_MANIFEST:-"${STAGED_DIR}/site_manifest.tsv"}
LOG_DIR=${LOG_DIR:-"${ANALYSIS_DIR}/logs/yellow_censoring_sites"}
ARRAY_INDEX=${FILE_INDEX:-${SLURM_ARRAY_TASK_ID:-1}}

AREA_RATIO=${AREA_RATIO:-2}
EXCLUDE_NEGATIVE_TRACK_IDS=${EXCLUDE_NEGATIVE_TRACK_IDS:-TRUE}
MIN_SEGMENT_FRAMES=${MIN_SEGMENT_FRAMES:-3}
OVERWRITE=${OVERWRITE:-FALSE}

mkdir -p "$STAGED_DIR" "$SITE_DIR" "$LOG_DIR"

if [ ! -s "$SITE_MANIFEST" ]; then
  MANIFEST_LOCK="${SITE_MANIFEST}.create.lock"
  if mkdir "$MANIFEST_LOCK" 2>/dev/null; then
    trap 'rm -rf "$MANIFEST_LOCK"' EXIT
    if [ ! -s "$SITE_MANIFEST" ]; then
      tmp_manifest="${SITE_MANIFEST}.tmp.$$"
      {
        printf "file_index\tsite_id\ttrackpoint_links\n"
        idx=0
        while IFS= read -r trackpoint_links; do
          base=$(basename "$trackpoint_links")
          site_id=${base%_nucleus_trackpoint_cpsam_links_framesall_excluding_frame0.tsv}
          idx=$((idx + 1))
          printf "%s\t%s\t%s\n" "$idx" "$site_id" "$trackpoint_links"
        done < <(find "$TRACKPOINT_LINKS_DIR" -maxdepth 1 -type f -name "*_nucleus_trackpoint_cpsam_links_framesall_excluding_frame0.tsv" | sort)
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

echo "[$(date --iso-8601=seconds)] Starting k00_yellow_censor task"
echo "HOSTNAME=${HOSTNAME:-unknown}"
echo "PWD=$(pwd)"
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS_DIR=${ANALYSIS_DIR}"
echo "TASK_SCRIPT=${TASK_SCRIPT}"
echo "TRACKPOINT_LINKS_DIR=${TRACKPOINT_LINKS_DIR}"
echo "NEAREST_DISTANCES_DIR=${NEAREST_DISTANCES_DIR}"
echo "SITE_MANIFEST=${SITE_MANIFEST}"
echo "SITE_DIR=${SITE_DIR}"
echo "MANIFEST_ROWS=${MANIFEST_ROWS}"
echo "ARRAY_INDEX=${ARRAY_INDEX}"
echo "SITE_ID=${SITE_ID:-}"
echo "AREA_RATIO=${AREA_RATIO}"
echo "MIN_SEGMENT_FRAMES=${MIN_SEGMENT_FRAMES}"
echo "EXCLUDE_NEGATIVE_TRACK_IDS=${EXCLUDE_NEGATIVE_TRACK_IDS}"
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
  --repo_root="${REPO_ROOT}" \
  --manifest="${SITE_MANIFEST}" \
  --file_index="${ARRAY_INDEX}" \
  --trackpoint_links_dir="${TRACKPOINT_LINKS_DIR}" \
  --nearest_distances_dir="${NEAREST_DISTANCES_DIR}" \
  --out_dir="${SITE_DIR}" \
  --area_ratio="${AREA_RATIO}" \
  --exclude_negative_track_ids="${EXCLUDE_NEGATIVE_TRACK_IDS}" \
  --min_segment_frames="${MIN_SEGMENT_FRAMES}" \
  --overwrite="${OVERWRITE}"

echo "[$(date --iso-8601=seconds)] Finished k00_yellow_censor task"
