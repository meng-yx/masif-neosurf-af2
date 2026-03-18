masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_root=$masif_neosurf_root/masif
masif_source=$masif_root/source/
docker_image=$masif_neosurf_root/masif-neosurf_v0.1.sif
export PYTHONPATH=$PYTHONPATH:$masif_source

SINGULARITY_BIND="$masif_neosurf_root:$masif_neosurf_root"

singularity exec --bind $SINGULARITY_BIND $docker_image python3 $masif_source/masif_ppi_search/masif_ppi_search_train.py $1