#!/bin/bash
#SBATCH --job-name=build_cache
#SBATCH --partition=l40s
#SBATCH --output=logs/build_cache-%A.out
#SBATCH --error=logs/build_cache-%A.out
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=32000


set -euo pipefail

usage() {
    cat <<'EOF'
Usage: sbatch build_cache.slurm [OPTIONS]

Options:
  --custom-params-module MODULE   Python module for custom params.
                                  Sets CUSTOM_PARAMS_MODULE.
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

export CUSTOM_PARAMS_MODULE=${CUSTOM_PARAMS_MODULE:-nn_models.v1.custom_params}

bash ./build_cache.sh "${CUSTOM_PARAMS_MODULE}"
