#!/bin/bash
set -euo pipefail

masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_root=$masif_neosurf_root/masif
masif_seed_root=$masif_neosurf_root/masif_seed_search
masif_source=$masif_root/source
masif_seed_source=$masif_seed_root/source
config_path=${1:-$masif_seed_root/data/protein_interaction_nn/configs/protein_interaction_v1.yaml}
split=${2:-test}
checkpoint=${3:-}

docker_image=$masif_neosurf_root/masif-neosurf_v0.1.sif
SINGULARITY_BIND="$masif_neosurf_root:$masif_neosurf_root"
export SINGULARITYENV_PYTHONPATH=$masif_source:$masif_seed_source
export SINGULARITYENV_LD_LIBRARY_PATH=/usr/local/lib/python3.6/site-packages/pymesh/lib:/usr/lib/x86_64-linux-gnu:/lib/x86_64-linux-gnu

extra_args=()
if [ -n "$checkpoint" ]; then
  extra_args+=(--checkpoint "$checkpoint")
fi

singularity exec --cleanenv --bind $SINGULARITY_BIND $docker_image \
  python3 -u "$masif_seed_source/protein_interaction/evaluate.py" \
  --config "$config_path" \
  --split "$split" \
  "${extra_args[@]}"

