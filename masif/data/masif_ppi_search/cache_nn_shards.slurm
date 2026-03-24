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

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: sbatch [--array=RANGE] cache_nn_shards.slurm [OPTIONS]

Options:
  --custom-params-module MODULE   Python module for custom params.
                                  Sets CUSTOM_PARAMS_MODULE.
  --catalog-dir DIR               Cache catalog directory.
                                  Sets CACHE_CATALOG_DIR.
  --subset-root DIR               Subset shard root directory.
                                  Sets CACHE_SUBSET_ROOT.
  --num-subsets INT               Number of subsets/shards.
                                  Sets NUM_SUBSETS.
  --pair-cache-max-records INT    Max records per in-memory pair cache.
                                  Sets PAIR_CACHE_MAX_RECORDS.
  --flush-every-pos INT           Flush frequency for positive samples.
                                  Sets SHARD_FLUSH_EVERY_POS.
  --subset-id INT                 Override subset id (arg1 to .sh). Useful
                                  for non-array debugging runs.
  -h, --help                      Show this help message.
EOF
}

subset_id_override=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --custom-params-module)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            export CUSTOM_PARAMS_MODULE="$2"
            shift 2
            ;;
        --catalog-dir)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            export CACHE_CATALOG_DIR="$2"
            shift 2
            ;;
        --subset-root)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            export CACHE_SUBSET_ROOT="$2"
            shift 2
            ;;
        --num-subsets)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            export NUM_SUBSETS="$2"
            shift 2
            ;;
        --pair-cache-max-records)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            export PAIR_CACHE_MAX_RECORDS="$2"
            shift 2
            ;;
        --flush-every-pos)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            export SHARD_FLUSH_EVERY_POS="$2"
            shift 2
            ;;
        --subset-id)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            subset_id_override="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

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
if [ -n "$subset_id_override" ]; then
    ./cache_nn_shards.sh "$subset_id_override"
else
    ./cache_nn_shards.sh
fi
