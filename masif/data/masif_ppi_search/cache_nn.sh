#!/bin/bash

# --- legacy monolithic cache generation ---
# Very slow and prone to OOM kill

masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_root=$masif_neosurf_root/masif
masif_source=$masif_root/source/
docker_image=$masif_neosurf_root/masif-neosurf_v0.1.sif
export PYTHONPATH=$PYTHONPATH:$masif_source

# Singularity bind paths - bind the entire repo root to make all paths accessible
SINGULARITY_BIND="$masif_neosurf_root:$masif_neosurf_root"

singularity exec --bind $SINGULARITY_BIND $docker_image python3 $masif_source/masif_ppi_search/masif_ppi_search_cache_training_data.py \
    nn_models.sc05.all_feat.custom_params