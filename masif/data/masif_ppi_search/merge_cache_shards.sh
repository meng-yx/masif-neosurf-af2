#!/bin/bash
set -euo pipefail

masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_root=$masif_neosurf_root/masif
masif_source=$masif_root/source/
docker_image=$masif_neosurf_root/masif-neosurf_v0.1.sif
export PYTHONPATH="${PYTHONPATH:-}:$masif_source"

SINGULARITY_BIND="$masif_neosurf_root:$masif_neosurf_root"
custom_params_module=${CUSTOM_PARAMS_MODULE:-nn_models.sc05.all_feat.custom_params}
subset_root=${CACHE_SUBSET_ROOT:-nn_models/sc05/cache/subsets}
cache_dir=${CACHE_OUTPUT_DIR:-nn_models/sc05/cache}
num_subsets=${NUM_SUBSETS:?NUM_SUBSETS must be set}

singularity exec --bind $SINGULARITY_BIND $docker_image \
    python3 $masif_source/masif_ppi_search/merge_cache_shards.py \
    $custom_params_module \
    --subset-root $subset_root \
    --num-subsets $num_subsets \
    --cache-dir $cache_dir
