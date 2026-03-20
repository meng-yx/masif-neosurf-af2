#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 7000
#SBATCH --time 01:00:00
#SBATCH --output=exelogs/build_cache_catalog-%A.out
#SBATCH --error=exelogs/build_cache_catalog-%A.out

./build_cache_catalog.sh
