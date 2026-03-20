#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 12G
#SBATCH --time 2:00:00
#SBATCH --partition=standard
#SBATCH --output=exelogs/eval_v1-%A.out
#SBATCH --error=exelogs/eval_v1-%A.out

echo "STARTING V1 EVALUATION AT $(date)"

CONFIG_PATH="${CONFIG_PATH:-configs/protein_interaction_v1.yaml}"
SPLIT="${SPLIT:-test}"
CHECKPOINT="${CHECKPOINT:-}"
./run_eval_v1.sh "$CONFIG_PATH" "$SPLIT" "$CHECKPOINT"

echo "FINISHED V1 EVALUATION AT $(date)"

