#!/bin/bash
#SBATCH --job-name=jupyter
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=logs/jupyter-%j.out
#SBATCH --error=logs/jupyter-%j.out

set -euo pipefail

echo "Starting Jupyter job on $(hostname)"
echo "SLURM_JOB_ID=$SLURM_JOB_ID"

# --- environment ---
IMAGE=/scratch/ymeng/tricomplex-design_Meng/masif-neosurf_v0.1.sif
REPO_ROOT=/scratch/ymeng/masif-neosurf-af2
SINGULARITY_BIND="$REPO_ROOT:$REPO_ROOT,/work:/work,/scratch:/scratch"

# --- thread control ---
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# --- pick a port ---
PORT=$(python -c "
import socket
s=socket.socket()
s.bind(('',0))
print(s.getsockname()[1])
s.close()
")

# --- launch Jupyter inside Singularity ---
singularity exec --cleanenv --bind $SINGULARITY_BIND $IMAGE \
  bash -c "pip install pandas numpy tqdm pyyaml && \
    jupyter notebook \
    --no-browser \
    --ip=0.0.0.0 \
    --port=${PORT} \
    --ServerApp.allow_remote_access=True \
    --notebook-dir=$REPO_ROOT"