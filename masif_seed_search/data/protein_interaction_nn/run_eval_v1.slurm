#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 12G
#SBATCH --time 2:00:00
#SBATCH --partition=standard
#SBATCH --output=exelogs/eval_v1-%A.out
#SBATCH --error=exelogs/eval_v1-%A.out

set -euo pipefail

submit_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"

echo "STARTING V1 EVALUATION AT $(date)"
bash "$submit_dir/run_eval_v1.sh" "$@"

echo "FINISHED V1 EVALUATION AT $(date)"

