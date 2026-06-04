#!/bin/bash
#SBATCH --job-name=masif_search
#SBATCH --output=logs/search-%A/slurm-%A_%a.out
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G

: '
# Submit from masif-neosurf repo root after writing numbered files under SEED_SUBSET_DIR:
#   sbatch --array=1-N scripts/slurm/search_array.sh <query_target> <out_dir> <seed_subset_dir>
#
# Example:
#   sbatch --array=1-5 scripts/slurm/search_array.sh 8VLB_A data/masif_search/8VLB_A data/masif_search/8VLB_A/subset
'

QUERY_TARGET=$1
OUT_DIR=$2
SEED_SUBSET_DIR=$3

if [[ -z "$QUERY_TARGET" || -z "$OUT_DIR" || -z "$SEED_SUBSET_DIR" ]]; then
    echo "Usage: sbatch --array=1-N $0 <query_target> <out_dir> <seed_subset_dir>" >&2
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

echo "Array task ${SLURM_ARRAY_TASK_ID}: ${QUERY_TARGET} subset ${SEED_SUBSET}"

bash scripts/bash/search.sh "${QUERY_TARGET}" "${OUT_DIR}" "${SEED_SUBSET}"
