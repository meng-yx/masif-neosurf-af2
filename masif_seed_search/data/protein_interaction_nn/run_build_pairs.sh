#!/bin/bash
set -euo pipefail

masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_seed_source=$masif_neosurf_root/masif_seed_search/source
config_path=${1:-$masif_neosurf_root/masif_seed_search/data/protein_interaction_nn/configs/protein_interaction_v1.yaml}
tiny_flag=${2:-""}

python_cmd="python3 -u $masif_seed_source/protein_interaction/build_pair_index.py --config $config_path"
if [ "$tiny_flag" = "--tiny" ]; then
  python_cmd="$python_cmd --tiny"
fi

PYTHONPATH=$masif_seed_source:$masif_neosurf_root/masif/source \
  bash -lc "$python_cmd"

