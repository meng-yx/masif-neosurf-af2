masif_root=../../../
masif_seed_root=../../
masif_source=$masif_root/source/
masif_seed_source=$masif_seed_root/source/
masif_matlab=$masif_root/source/matlab_libs/
masif_data=$masif_root/data/
export PYTHONPATH=$PYTHONPATH:$masif_source:$masif_seed_source

# Construct absolute path to PDB list file
# Note: Change to testing_seed_benchmark.txt if generating testing data instead of training data
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIST_FILE="${SCRIPT_DIR}/lists/training_seed_benchmark.txt"

python -W ignore -u $masif_seed_source/seed_search_generate_training_data.py \
  --data-dir $masif_root/data/masif_ppi_search/ \
  --K 1000 \
  --ransac-iter 2000 \
  --patch-radius 9 \
  --output-dir training_data_9A/ \
  --pdb-list-index $1 \
  --pdb-list "$LIST_FILE"
