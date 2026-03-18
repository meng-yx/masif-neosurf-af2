#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 12G
#SBATCH --time 01:00:00
#SBATCH --output=exelogs/cache_nn.%A_%a.out
#SBATCH --error=exelogs/cache_nn.%A_%a.out

./cache_nn.sh