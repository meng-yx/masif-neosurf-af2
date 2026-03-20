#!/bin/bash
set -euo pipefail

masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_root=$masif_neosurf_root/masif
masif_seed_root=$masif_neosurf_root/masif_seed_search
masif_source=$masif_root/source
masif_seed_source=$masif_seed_root/source
default_config="$masif_seed_root/data/protein_interaction_nn/configs/protein_interaction_v1.yaml"
config_path="$default_config"
split="test"
checkpoint=""
pairs_csv=""

usage() {
  cat <<'EOF'
Usage:
  run_eval_v1.sh [--config <path>] [--split <name>] [--checkpoint <path>] [--pairs-csv <path>]

Options:
  --config <path>      Config YAML/JSON path.
  --split <name>       Split to evaluate: train|val|test (default: test).
  --checkpoint <path>  Optional checkpoint override.
  --pairs-csv <path>   Optional pairs CSV override.
  -h, --help           Show this help message.
EOF
}

# Backward compatibility for old positional usage:
#   run_eval_v1.sh <config_path> [split] [checkpoint]
if [ "${1:-}" != "" ] && [[ "${1:-}" != --* ]]; then
  config_path="$1"
  shift
  if [ "${1:-}" != "" ] && [[ "${1:-}" != --* ]]; then
    split="$1"
    shift
  fi
  if [ "${1:-}" != "" ] && [[ "${1:-}" != --* ]]; then
    checkpoint="$1"
    shift
  fi
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
    --split)
      if [[ $# -lt 2 ]]; then
        echo "Error: --split requires a value" >&2
        usage
        exit 2
      fi
      split="$2"
      shift 2
      ;;
    --checkpoint)
      if [[ $# -lt 2 ]]; then
        echo "Error: --checkpoint requires a value" >&2
        usage
        exit 2
      fi
      checkpoint="$2"
      shift 2
      ;;
    --pairs-csv)
      if [[ $# -lt 2 ]]; then
        echo "Error: --pairs-csv requires a value" >&2
        usage
        exit 2
      fi
      pairs_csv="$2"
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
if [ -n "$checkpoint" ]; then
  checkpoint=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$checkpoint")
fi
if [ -n "$pairs_csv" ]; then
  pairs_csv=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$pairs_csv")
fi

docker_image=$masif_neosurf_root/masif-neosurf_v0.1.sif
SINGULARITY_BIND="$masif_neosurf_root:$masif_neosurf_root"
export SINGULARITYENV_PYTHONPATH=$masif_source:$masif_seed_source
export SINGULARITYENV_LD_LIBRARY_PATH=/usr/local/lib/python3.6/site-packages/pymesh/lib:/usr/lib/x86_64-linux-gnu:/lib/x86_64-linux-gnu

extra_args=()
if [ -n "$checkpoint" ]; then
  extra_args+=(--checkpoint "$checkpoint")
fi
if [ -n "$pairs_csv" ]; then
  extra_args+=(--pairs_csv "$pairs_csv")
fi

singularity exec --cleanenv --bind $SINGULARITY_BIND $docker_image \
  python3 -u "$masif_seed_source/protein_interaction/evaluate.py" \
  --config "$config_path" \
  --split "$split" \
  "${extra_args[@]}"

