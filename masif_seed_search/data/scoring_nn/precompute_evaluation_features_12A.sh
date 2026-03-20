#!/bin/bash

# Legacy code without singularity container
# masif_neosurf_root=$(git rev-parse --show-toplevel)
# masif_root=$masif_neosurf_root/masif
# masif_source=$masif_root/source/
# masif_seed_source=../../source/
# masif_data=$masif_root/masif/data/
# export PYTHONPATH=$PYTHONPATH:$masif_source:
# python -W ignore $masif_seed_source/precompute_evaluation_features.py training_data_12A_seed_benchmark/ $1


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

# Construct absolute path to PDB list file
# Note: Change to testing_seed_benchmark.txt if generating testing data instead of training data
LIST_FILE="${masif_seed_root}/data/scoring_nn/lists/full_list.txt"
DATA_DIR="${SCORING_FEATURES_INPUT_DIR:-training_data_12A_seed_benchmark/}"

if [ ! -d "${DATA_DIR}$1" ]; then
  echo "ERROR: Expected input directory not found: ${DATA_DIR}$1"
  echo "Hint: ensure make_transformations output dir matches precompute input dir."
  exit 2
fi

singularity exec --cleanenv --bind $SINGULARITY_BIND $docker_image python -W ignore -u $masif_seed_source/precompute_evaluation_features.py "$DATA_DIR" $1