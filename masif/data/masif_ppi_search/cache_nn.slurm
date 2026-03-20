#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=bigmem
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 64
#SBATCH --mem 896000
#SBATCH --time 08:00:00
#SBATCH --output=exelogs/cache_nn.%A_%a.out
#SBATCH --error=exelogs/cache_nn.%A_%a.out

./cache_nn.sh