#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 7000
#SBATCH --time 08:00:00
#SBATCH --output=exelogs/cache_nn.%A_%a.out
#SBATCH --error=exelogs/cache_nn.%A_%a.out

# --- legacy monolithic cache generation ---
# Very slow and prone to OOM kill

./cache_nn.sh