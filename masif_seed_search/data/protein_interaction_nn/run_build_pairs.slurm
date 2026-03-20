#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 8G
#SBATCH --time 1:00:00
#SBATCH --partition=standard
#SBATCH --output=exelogs/build_pairs-%A.out
#SBATCH --error=exelogs/build_pairs-%A.out

set -euo pipefail

submit_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"

echo "STARTING PAIR INDEX BUILD AT $(date)"
bash "$submit_dir/run_build_pairs.sh" "$@"
echo "FINISHED PAIR INDEX BUILD AT $(date)"
