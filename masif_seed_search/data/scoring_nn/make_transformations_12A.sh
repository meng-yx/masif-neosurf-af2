#!/bin/bash
masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_root=$masif_neosurf_root/masif
masif_seed_root=$masif_neosurf_root/masif_seed_search
masif_source=$masif_root/source/
masif_seed_source=$masif_seed_root/source

# Run inside singularity container
docker_image=$masif_neosurf_root/masif-neosurf_v0.1.sif
SINGULARITY_BIND="$masif_neosurf_root:$masif_neosurf_root"
# Set PYTHONPATH inside container (SINGULARITYENV_ prefix passes env vars into container)
export SINGULARITYENV_PYTHONPATH=$masif_source:$masif_seed_source
# Set LD_LIBRARY_PATH to prioritize pymesh's bundled libstdc++ which has GLIBCXX_3.4.26
# The system libstdc++ only has up to GLIBCXX_3.4.25, causing the mismatch error
# Prepend pymesh lib directory so its libstdc++.so.6 is found first
export SINGULARITYENV_LD_LIBRARY_PATH=/usr/local/lib/python3.6/site-packages/pymesh/lib:/usr/lib/x86_64-linux-gnu:/lib/x86_64-linux-gnu
# Use --cleanenv to prevent host environment from interfering with container libraries
singularity exec --cleanenv --bind $SINGULARITY_BIND $docker_image python -W ignore -u $masif_seed_source/seed_search_generate_training_data.py $masif_root/data/masif_ppi_search/ 1000 2000 12 testing_data_12A_seed_benchmark/ $1
