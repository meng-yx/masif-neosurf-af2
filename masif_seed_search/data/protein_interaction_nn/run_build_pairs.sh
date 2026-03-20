#!/bin/bash
set -euo pipefail

masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_seed_source=$masif_neosurf_root/masif_seed_search/source
default_config="$masif_neosurf_root/masif_seed_search/data/protein_interaction_nn/configs/protein_interaction_v1.yaml"
config_path="$default_config"
tiny_enabled=0

usage() {
  cat <<'EOF'
Usage:
  run_build_pairs.sh [--config <path>] [--tiny]

Options:
  --config <path>  Config YAML/JSON path.
  --tiny           Enable tiny subset mode.
  -h, --help       Show this help message.
EOF
}

# Backward compatibility for old positional usage:
#   run_build_pairs.sh <config_path> [--tiny]
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
    --tiny)
      tiny_enabled=1
      shift
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

args=(
  python3
  -u
  "$masif_seed_source/protein_interaction/build_pair_index.py"
  --config
  "$config_path"
)
if [[ "$tiny_enabled" -eq 1 ]]; then
  args+=(--tiny)
fi

PYTHONPATH=$masif_seed_source:$masif_neosurf_root/masif/source \
  "${args[@]}"

