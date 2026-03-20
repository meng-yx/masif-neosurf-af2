#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 24G
#SBATCH --time 8:00:00
#SBATCH --partition=standard
#SBATCH --output=exelogs/train_v1-%A.out
#SBATCH --error=exelogs/train_v1-%A.out

set -euo pipefail

submit_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"

echo "STARTING V1 TRAINING AT $(date)"
bash "$submit_dir/run_train_v1.sh" "$@"

echo "FINISHED V1 TRAINING AT $(date)"

