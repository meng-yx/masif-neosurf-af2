#!/bin/bash
#SBATCH --job-name=masif_search
#SBATCH --output=logs/search-%A/slurm-%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7000

: '
# Submit from masif-neosurf repo root after writing numbered files under SEED_SUBSET_DIR.
#
# Single-target mode (original):
#   sbatch --array=1-N scripts/slurm/search_array.sh <query_target> <out_dir> <seed_subset_dir>
#
# Example (single-target):
#   sbatch --array=1-5 scripts/slurm/search_array.sh 8VLB_A data/masif_search data/masif_search/subset
#
# Multi-target mode (pass a .txt manifest as the first argument):
#   sbatch --array=1-N scripts/slurm/search_array.sh <query_targets.txt> <masif_search_out_dir> <seed_subset_dir>
#
# Each array task reads the shared seed subset file and loops over every target
# in query_targets.txt, writing results to <masif_search_out_dir>/<target>/.
#
# Example (multi-target):
#   sbatch --array=1-500 scripts/slurm/search_array.sh \
#       data/masif_search/query_targets.txt \
#       data/masif_search \
#       data/masif_search/subset
'

QUERY_TARGET=$1
OUT_DIR=$2
SEED_SUBSET_DIR=$3

if [[ -z "$QUERY_TARGET" || -z "$OUT_DIR" || -z "$SEED_SUBSET_DIR" ]]; then
    echo "Usage: sbatch --array=1-N $0 <query_target_or_txt> <out_dir> <seed_subset_dir>" >&2
    exit 1
fi

if [[ ! -f scripts/config.sh ]]; then
    echo "Error: run sbatch from masif-neosurf repo root (scripts/config.sh not found in $(pwd))" >&2
    exit 1
fi

source scripts/config.sh

: "${SLURM_ARRAY_TASK_ID:=1}"
SEED_SUBSET="${SEED_SUBSET_DIR}/${SLURM_ARRAY_TASK_ID}"

if [[ ! -f "${SEED_SUBSET}" ]]; then
    echo "Error: seed subset file not found: ${SEED_SUBSET}" >&2
    exit 1
fi

if [[ -f "${QUERY_TARGET}" && "${QUERY_TARGET}" == *.txt ]]; then
    # Multi-target mode: loop over every target in the manifest file
    while IFS= read -r target || [[ -n "${target}" ]]; do
        [[ -z "${target}" || "${target}" =~ ^# ]] && continue
        echo "Array task ${SLURM_ARRAY_TASK_ID}: ${target} subset ${SEED_SUBSET}"
        if ! bash scripts/bash/search.sh "${target}" "${OUT_DIR}" "${SEED_SUBSET}"; then
            echo "Error: search failed for ${target}" >&2
        fi
    done < "${QUERY_TARGET}"
else
    # Single-target mode (original behaviour)
    echo "Array task ${SLURM_ARRAY_TASK_ID}: ${QUERY_TARGET} subset ${SEED_SUBSET}"
    bash scripts/bash/search.sh "${QUERY_TARGET}" "${OUT_DIR}" "${SEED_SUBSET}"
fi
