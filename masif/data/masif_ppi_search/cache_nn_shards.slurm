#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 7000
#SBATCH --time 01:00:00
#SBATCH --output=exelogs/cache_nn_shards-%A/slurm-%A_%a.out
#SBATCH --error=exelogs/cache_nn_shards-%A/slurm-%A_%a.out

# submit array job as:
# sbatch --array=0-999 cache_nn_shards.slurm

# Check of subset _SUCCESS file already exists
# - Use this to skip already completed subsets to avoid re-running the same subset
# if [ -f "nn_models/sc05/cache/subsets/subset_$SLURM_ARRAY_TASK_ID/_SUCCESS" ]; then
#     echo "Subset $SLURM_ARRAY_TASK_ID already exists"
#     exit 0
# fi

# Keep this in sync with the --array range.
export NUM_SUBSETS=${NUM_SUBSETS:-1000}
export PAIR_CACHE_MAX_RECORDS=${PAIR_CACHE_MAX_RECORDS:-16}
export SHARD_FLUSH_EVERY_POS=${SHARD_FLUSH_EVERY_POS:-256}
echo "NUM_SUBSETS: $NUM_SUBSETS"
echo "PAIR_CACHE_MAX_RECORDS: $PAIR_CACHE_MAX_RECORDS"
echo "SHARD_FLUSH_EVERY_POS: $SHARD_FLUSH_EVERY_POS"
./cache_nn_shards.sh
