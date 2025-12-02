#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 16G
#SBATCH --time 4:00:00
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu

echo "STARTING TRAINING AT $(date)"

# Load your environment
# source /path/to/your/environment.sh

masif_root=../../../
masif_seed_root=../../
masif_source=$masif_root/masif/source/
masif_seed_source=$masif_seed_root/source/
export PYTHONPATH=$PYTHONPATH:$masif_source:$masif_seed_source

# Train with all 4 features for 12A patches
python3 $masif_seed_source/train_alignment_evaluation_nn.py 12A 0_1_2_3

echo "FINISHED TRAINING AT $(date)"