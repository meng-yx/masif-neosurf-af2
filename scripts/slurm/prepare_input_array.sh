#!/bin/bash
#SBATCH --job-name=prepare_input
#SBATCH --output=logs/prepare_input-%A/slurm-%A_%a.out
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7000

: '
# Submit from repo root after writing numbered subset CSVs under INPUT_SUBSET_DIR:
#   sbatch --array=1-N scripts/slurm/prepare_input_array.sh \
#       <input_subset_dir> <outdir> <output_subset_dir> <evoef2_bin>
#
# Array task k reads  input_subsets/input_k.csv
#              writes output_subsets/output_k.csv
#              writes .pdb/.sdf files to outdir (shared across tasks)
#
# Example:
#   sbatch --array=1-100 scripts/slurm/prepare_input_array.sh \
#       data/processing/1_prep_input/input_subsets \
#       data/input \
#       data/processing/1_prep_input/output_subsets \
#       EvoEF2/EvoEF2
'

INPUT_SUBSET_DIR=$1
OUTDIR=$2
OUTPUT_SUBSET_DIR=$3
EVOEF2_BIN=$4

if [[ -z "$INPUT_SUBSET_DIR" || -z "$OUTDIR" || -z "$OUTPUT_SUBSET_DIR" || -z "$EVOEF2_BIN" ]]; then
    echo "Usage: sbatch --array=1-N $0 <input_subset_dir> <outdir> <output_subset_dir> <evoef2_bin>" >&2
    exit 1
fi

if [[ ! -f scripts/config.sh ]]; then
    echo "Error: run sbatch from repo root (scripts/config.sh not found in $(pwd))" >&2
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

if [[ ! -f "${EVOEF2_BIN}" ]]; then
    echo "Error: EvoEF2 binary not found: ${EVOEF2_BIN}" >&2
    exit 1
fi

mkdir -p "${OUTDIR}" "${OUTPUT_SUBSET_DIR}" "logs/prepare_input-${SLURM_ARRAY_JOB_ID:-local}"

echo "Array task ${SLURM_ARRAY_TASK_ID}: ${INPUT_CSV} -> ${OUTPUT_CSV}"

source ~/.bashrc
conda activate MaSIF

python scripts/python/prepare_input.py \
    --input_csv "${INPUT_CSV}" \
    --outdir "${OUTDIR}" \
    --out_csv "${OUTPUT_CSV}" \
    --evoef2_bin "${EVOEF2_BIN}"
