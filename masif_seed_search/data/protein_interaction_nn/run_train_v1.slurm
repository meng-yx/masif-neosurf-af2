#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 24G
#SBATCH --time 8:00:00
#SBATCH --partition=standard
#SBATCH --output=exelogs/train_v1-%A.out
#SBATCH --error=exelogs/train_v1-%A.out

echo "STARTING V1 TRAINING AT $(date)"

CONFIG_PATH="${CONFIG_PATH:-configs/protein_interaction_v1.yaml}"
./run_train_v1.sh "$CONFIG_PATH"

echo "FINISHED V1 TRAINING AT $(date)"

