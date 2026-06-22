#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --output=logs/preprocess-%A/slurm-%A_%a.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=14000

: '
# Submit from masif-neosurf repo root after writing numbered subset CSVs:
#   sbatch --array=1-N scripts/slurm/preprocess_array.sh \
#       <input_subset_dir> <output_subset_dir>
#
# Array task k reads  input_subsets/input_k.csv
#              writes output_subsets/output_k.csv
#
# Example:
#   sbatch --array=1-100 scripts/slurm/preprocess_array.sh \
#       data/processing/2_masif_preprocess/input_subsets \
#       data/processing/2_masif_preprocess/output_subsets
'

INPUT_SUBSET_DIR=$1
OUTPUT_SUBSET_DIR=$2

if [[ -z "$INPUT_SUBSET_DIR" || -z "$OUTPUT_SUBSET_DIR" ]]; then
    echo "Usage: sbatch --array=1-N $0 <input_subset_dir> <output_subset_dir>" >&2
    exit 1
fi

if [[ ! -f scripts/config.sh ]]; then
    echo "Error: run sbatch from masif-neosurf repo root (scripts/config.sh not found in $(pwd))" >&2
    exit 1
fi

source scripts/config.sh

: "${SLURM_ARRAY_TASK_ID:=1}"

INPUT_CSV="${INPUT_SUBSET_DIR}/input_${SLURM_ARRAY_TASK_ID}.csv"
OUTPUT_CSV="${OUTPUT_SUBSET_DIR}/output_${SLURM_ARRAY_TASK_ID}.csv"

if [[ ! -f "${INPUT_CSV}" ]]; then
    echo "No input subset for array task ${SLURM_ARRAY_TASK_ID}: ${INPUT_CSV}; skipping."
    exit 0
fi

mkdir -p "${OUTPUT_SUBSET_DIR}" "logs/preprocess-${SLURM_ARRAY_JOB_ID:-local}"

echo "Array task ${SLURM_ARRAY_TASK_ID}: ${INPUT_CSV} -> ${OUTPUT_CSV}"

source ~/.bashrc
conda activate MaSIF

PREPROCESS_ARGS=(
    --input_csv "${INPUT_CSV}"
    --out_csv "${OUTPUT_CSV}"
    --preprocess_root "${TARGET_PREPROCESS_ROOT}"
    --repo_root "$(pwd)"
)

python scripts/python/preprocess_manifest.py "${PREPROCESS_ARGS[@]}"
