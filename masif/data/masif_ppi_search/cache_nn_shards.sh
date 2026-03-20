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
subset_root=${CACHE_SUBSET_ROOT:-nn_models/sc05/cache/subsets}
num_subsets=${NUM_SUBSETS:?NUM_SUBSETS must be set}
subset_id=${1:-${SLURM_ARRAY_TASK_ID:-}}
pair_cache_max_records=${PAIR_CACHE_MAX_RECORDS:-16}
flush_every_pos=${SHARD_FLUSH_EVERY_POS:-256}

if [ -z "${subset_id}" ]; then
    echo "subset id is not provided (arg1 or SLURM_ARRAY_TASK_ID)" >&2
    exit 1
fi

singularity exec --bind $SINGULARITY_BIND $docker_image \
    python3 $masif_source/masif_ppi_search/masif_ppi_search_cache_training_data.py \
    $custom_params_module \
    --catalog-dir $catalog_dir \
    --subset-root $subset_root \
    --subset-id $subset_id \
    --num-subsets $num_subsets \
    --pair-cache-max-records $pair_cache_max_records \
    --flush-every-pos $flush_every_pos
