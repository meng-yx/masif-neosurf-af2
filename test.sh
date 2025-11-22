#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 7000
#SBATCH --time 01:00:00
#SBATCH --array=1
#SBATCH --output=exelogs/_masif_precompute.%A_%a.out
#SBATCH --error=exelogs/_masif_precompute.%A_%a.out

echo 'Hello World'