#!/bin/bash
set -euo pipefail

masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_root=$masif_neosurf_root/masif
masif_seed_root=$masif_neosurf_root/masif_seed_search
masif_source=$masif_root/source
masif_seed_source=$masif_seed_root/source
default_config="$masif_seed_root/data/protein_interaction_nn/configs/protein_interaction_v1.yaml"
config_path="$default_config"
train_pairs=""
val_pairs=""

usage() {
  cat <<'EOF'
Usage:
  run_train_v1.sh [--config <path>] [--train-pairs <csv>] [--val-pairs <csv>]

Options:
  --config <path>      Config YAML/JSON path.
  --train-pairs <csv>  Optional train pairs CSV override.
  --val-pairs <csv>    Optional val pairs CSV override.
  -h, --help           Show this help message.
EOF
}

# Backward compatibility for old positional usage:
#   run_train_v1.sh <config_path>
if [ "${1:-}" != "" ] && [[ "${1:-}" != --* ]]; then
  config_path="$1"
  shift
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      if [[ $# -lt 2 ]]; then
        echo "Error: --config requires a value" >&2
        usage
        exit 2
      fi
      config_path="$2"
      shift 2
      ;;
    --train-pairs)
      if [[ $# -lt 2 ]]; then
        echo "Error: --train-pairs requires a value" >&2
        usage
        exit 2
      fi
      train_pairs="$2"
      shift 2
      ;;
    --val-pairs)
      if [[ $# -lt 2 ]]; then
        echo "Error: --val-pairs requires a value" >&2
        usage
        exit 2
      fi
      val_pairs="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

config_path=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$config_path")
if [ -n "$train_pairs" ]; then
  train_pairs=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$train_pairs")
fi
if [ -n "$val_pairs" ]; then
  val_pairs=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$val_pairs")
fi

docker_image=$masif_neosurf_root/masif-neosurf_v0.1.sif
SINGULARITY_BIND="$masif_neosurf_root:$masif_neosurf_root"
export SINGULARITYENV_PYTHONPATH=$masif_source:$masif_seed_source
export SINGULARITYENV_LD_LIBRARY_PATH=/usr/local/lib/python3.6/site-packages/pymesh/lib:/usr/lib/x86_64-linux-gnu:/lib/x86_64-linux-gnu

extra_args=()
if [ -n "$train_pairs" ]; then
  extra_args+=(--train_pairs "$train_pairs")
fi
if [ -n "$val_pairs" ]; then
  extra_args+=(--val_pairs "$val_pairs")
fi

singularity exec --cleanenv --bind $SINGULARITY_BIND $docker_image \
  python3 -u "$masif_seed_source/protein_interaction/train.py" \
  --config "$config_path" \
  "${extra_args[@]}"

