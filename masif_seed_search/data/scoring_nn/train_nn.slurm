#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 16G
#SBATCH --time 4:00:00
#SBATCH --partition=standard
#SBATCH --output=exelogs/train_nn-%A.out
#SBATCH --error=exelogs/train_nn-%A.out

echo "STARTING TRAINING AT $(date)"

# Train with defaults from train_nn.sh unless overridden.
PATCH_SIZE="${PATCH_SIZE:-12A}"
SELECTED_FEATURES="${SELECTED_FEATURES:-0_1_2_3}"
./train_nn.sh "$PATCH_SIZE" "$SELECTED_FEATURES"

echo "FINISHED TRAINING AT $(date)"