#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 28000
#SBATCH --time 01:00:00
#SBATCH --output=exelogs/merge_cache_shards-%A.out
#SBATCH --error=exelogs/merge_cache_shards-%A.out

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: sbatch merge_cache_shards.slurm [OPTIONS]

Options:
  --custom-params-module MODULE   Python module for custom params.
                                  Sets CUSTOM_PARAMS_MODULE.
  --subset-root DIR               Subset shard root directory.
                                  Sets CACHE_SUBSET_ROOT.
  --cache-output-dir DIR          Final merged cache output directory.
                                  Sets CACHE_OUTPUT_DIR.
  --num-subsets INT               Number of subsets/shards to merge.
                                  Sets NUM_SUBSETS.
  -h, --help                      Show this help message.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --custom-params-module)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            export CUSTOM_PARAMS_MODULE="$2"
            shift 2
            ;;
        --subset-root)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            export CACHE_SUBSET_ROOT="$2"
            shift 2
            ;;
        --cache-output-dir)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            export CACHE_OUTPUT_DIR="$2"
            shift 2
            ;;
        --num-subsets)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            export NUM_SUBSETS="$2"
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

export NUM_SUBSETS=${NUM_SUBSETS:-1000}
./merge_cache_shards.sh
