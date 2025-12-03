#!/bin/bash
#SBATCH --account=lpdi
#SBATCH --nodes 1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 7000
#SBATCH --time 01:00:00
#SBATCH --array=1
#SBATCH --output=exelogs/cache_nn.%A_%a.out
#SBATCH --error=exelogs/cache_nn.%A_%a.out

./cache_nn.sh
