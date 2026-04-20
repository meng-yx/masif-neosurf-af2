#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=00:24:00
#SBATCH --partition=h100
#SBATCH --qos=long
#SBATCH --gpus-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/train_nn.%A.out
#SBATCH --error=logs/train_nn.%A.out

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: sbatch train_nn.slurm [OPTIONS]

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


./train_nn.sh "${CUSTOM_PARAMS_MODULE}"
