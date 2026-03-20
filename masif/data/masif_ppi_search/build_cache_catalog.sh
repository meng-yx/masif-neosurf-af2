#!/bin/bash
set -euo pipefail

masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_root=$masif_neosurf_root/masif
masif_source=$masif_root/source/
docker_image=$masif_neosurf_root/masif-neosurf_v0.1.sif
export PYTHONPATH="${PYTHONPATH:-}:$masif_source"

SINGULARITY_BIND="$masif_neosurf_root:$masif_neosurf_root"
custom_params_module=${CUSTOM_PARAMS_MODULE:-nn_models.sc05.all_feat.custom_params}
catalog_dir=${CACHE_CATALOG_DIR:-nn_models/sc05/cache/catalog}
seed_arg=""
if [ -n "${CACHE_CATALOG_SEED:-}" ]; then
    seed_arg="--seed ${CACHE_CATALOG_SEED}"
fi

singularity exec --bind $SINGULARITY_BIND $docker_image \
    python3 $masif_source/masif_ppi_search/build_cache_catalog.py \
    $custom_params_module \
    --catalog-dir $catalog_dir \
    $seed_arg
