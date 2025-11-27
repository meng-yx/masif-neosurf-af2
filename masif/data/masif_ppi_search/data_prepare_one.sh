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

# Run Python scripts inside Singularity container from the masif/source directory
# This ensures relative paths in masif_opts resolve correctly
#singularity exec --bind $SINGULARITY_BIND $docker_image python 00a-pdb_to_af2.py $1

#singularity exec --bind $SINGULARITY_BIND $docker_image python $masif_source/data_preparation/00-pdb_download.py $1 
singularity exec --bind $SINGULARITY_BIND $docker_image python $masif_source/data_preparation/01-pdb_extract_and_triangulate.py $PDB_ID\_$CHAIN1 
singularity exec --bind $SINGULARITY_BIND $docker_image python $masif_source/data_preparation/01-pdb_extract_and_triangulate.py $PDB_ID\_$CHAIN2
singularity exec --bind $SINGULARITY_BIND $docker_image python $masif_source/data_preparation/04-masif_precompute.py masif_site $1
singularity exec --bind $SINGULARITY_BIND $docker_image python $masif_source/data_preparation/04-masif_precompute.py masif_ppi_search $1
