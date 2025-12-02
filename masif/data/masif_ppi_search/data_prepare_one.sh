masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_root=$masif_neosurf_root/masif
masif_source=$masif_root/source/
docker_image=$masif_neosurf_root/masif-neosurf_v0.1.sif
export PYTHONPATH=$PYTHONPATH:$masif_source
PDB_ID=$(echo $1| cut -d"_" -f1)
CHAIN1=$(echo $1| cut -d"_" -f2)
CHAIN2=$(echo $1| cut -d"_" -f3)

# Singularity bind paths - bind the entire repo root to make all paths accessible
SINGULARITY_BIND="$masif_neosurf_root:$masif_neosurf_root"

# The new AF2 download script uses conda environment
eval "$(conda shell.bash hook)"
conda activate MaSIF
python $masif_source/data_preparation/00a-PDB_to_AF2.py $1 \
     --pdb-out-dir ./data_preparation/00-raw_pdbs/ \
     --rmsd-out-dir ./data_preparation/00a-rmsd/ \
     --temp-dir ./data_preparation/temp/

# Run Python scripts inside Singularity container from the masif/source directory
# This ensures relative paths in masif_opts resolve correctly
singularity exec --bind $SINGULARITY_BIND $docker_image python $masif_source/data_preparation/01-pdb_extract_and_triangulate.py $PDB_ID\_$CHAIN1
singularity exec --bind $SINGULARITY_BIND $docker_image python $masif_source/data_preparation/01-pdb_extract_and_triangulate.py $PDB_ID\_$CHAIN2
singularity exec --bind $SINGULARITY_BIND $docker_image python $masif_source/data_preparation/04-masif_precompute.py masif_site $1
singularity exec --bind $SINGULARITY_BIND $docker_image python $masif_source/data_preparation/04-masif_precompute.py masif_ppi_search $1
