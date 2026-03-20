#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 28000
#SBATCH --time 01:00:00
#SBATCH --output=exelogs/merge_cache_shards-%A.out
#SBATCH --error=exelogs/merge_cache_shards-%A.out

export NUM_SUBSETS=${NUM_SUBSETS:-100}
./merge_cache_shards.sh
