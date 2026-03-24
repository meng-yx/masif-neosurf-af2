#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 7000
#SBATCH --time 01:00:00
#SBATCH --output=exelogs/build_cache_catalog-%A.out
#SBATCH --error=exelogs/build_cache_catalog-%A.out

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: sbatch build_cache_catalog.slurm [OPTIONS]

Options:
  --custom-params-module MODULE   Python module for custom params.
                                  Sets CUSTOM_PARAMS_MODULE.
  --catalog-dir DIR               Output catalog directory.
                                  Sets CACHE_CATALOG_DIR.
  --seed INT                      Optional random seed.
                                  Sets CACHE_CATALOG_SEED.
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
        --catalog-dir)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            export CACHE_CATALOG_DIR="$2"
            shift 2
            ;;
        --seed)
            [[ $# -ge 2 ]] || { echo "Missing value for $1" >&2; exit 2; }
            export CACHE_CATALOG_SEED="$2"
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

./build_cache_catalog.sh
