masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_root=$masif_neosurf_root/masif
masif_source=$masif_root/source/
docker_image=$masif_neosurf_root/masif-neosurf_v0.1.sif
echo "docker image: $docker_image"
export PYTHONPATH=$PYTHONPATH:$masif_source
if [ -z "$2" ]; then
    echo "Usage: $0 PDBID_CHAIN1_CHAIN2 /path/to/pinder_pdb_dir"
    exit 1
fi
MODEL=$(echo $1| cut -d"_" -f1)
CHAIN1=$(echo $1| cut -d"_" -f2)
CHAIN2=$(echo $1| cut -d"_" -f3)
PDB_ID=$(echo $MODEL | cut -d"-" -f1)
MODEL_TYPE=$(echo $1 | cut -d"-" -f2 | cut -d"_" -f1)
PINDER_PDB_DIR=$2
PINDER_PDB_DIR_REAL=$(readlink -f "$PINDER_PDB_DIR")

if [ ! -d "$PINDER_PDB_DIR_REAL" ]; then
    echo "Error: pinder_pdb_dir does not exist or is not a directory: $PINDER_PDB_DIR"
    exit 1
fi

# Bind repo root and the resolved Pinder directory (important when pinderMaSIF is a symlink).
SINGULARITY_BIND="$masif_neosurf_root:$masif_neosurf_root,$PINDER_PDB_DIR_REAL:$PINDER_PDB_DIR_REAL"

# Construct the path to the expected final output .npy file
FINAL_OUT_PATH_1="data_preparation/04b-precomputation_12A/precomputation/${MODEL}_${CHAIN1}_${CHAIN2}/p1_sc_labels.npy"
FINAL_OUT_PATH_2="data_preparation/04b-precomputation_12A/precomputation/${MODEL}_${CHAIN1}_${CHAIN2}/p2_sc_labels.npy"
if [ -f "$FINAL_OUT_PATH_1" ] && [ -f "$FINAL_OUT_PATH_2" ]; then
    echo "$FINAL_OUT_PATH_1 and $FINAL_OUT_PATH_2 already exist"
    echo "Skipped data preparation"
    exit 0
fi

echo "Running 00d-pdb_from_pinder.py"
# Run Python scripts inside Singularity container
echo "Using Pinder PDB directory: $PINDER_PDB_DIR_REAL"


# Pinder-specific logic for handling aligned AF2 models (AP/PA)
# ---- New: check if the complex_id is a PP model -----
if [ "$MODEL_TYPE" == "PP" ]; then
    echo "PP model detected, only prepare the PP model"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/00d-pdb_from_pinder.py" "$1" --pinder-pdb-dir "$PINDER_PDB_DIR_REAL"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/01-pdb_extract_and_triangulate.py" "${MODEL}_${CHAIN1}"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/01-pdb_extract_and_triangulate.py" "${MODEL}_${CHAIN2}"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/04-masif_precompute.py" masif_site "$1"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/04-masif_precompute.py" masif_ppi_search "$1"
    exit 0
else
    echo "Non-PP model detected, prepare the PP model first..."
    echo "Equivalent PP model: ${PDB_ID}-PP_${CHAIN1}_${CHAIN2}"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/00d-pdb_from_pinder.py" "${PDB_ID}-PP_${CHAIN1}_${CHAIN2}" --pinder-pdb-dir "$PINDER_PDB_DIR_REAL"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/01-pdb_extract_and_triangulate.py" "${PDB_ID}-PP_${CHAIN1}"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/01-pdb_extract_and_triangulate.py" "${PDB_ID}-PP_${CHAIN2}"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/04-masif_precompute.py" masif_site "${PDB_ID}-PP_${CHAIN1}_${CHAIN2}"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/04-masif_precompute.py" masif_ppi_search "${PDB_ID}-PP_${CHAIN1}_${CHAIN2}"
    echo "Preparing ${1}..."
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/00d-pdb_from_pinder.py" "$1" --pinder-pdb-dir "$PINDER_PDB_DIR_REAL"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/01-pdb_extract_and_triangulate.py" "${MODEL}_${CHAIN1}"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/01-pdb_extract_and_triangulate.py" "${MODEL}_${CHAIN2}"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/04-masif_precompute.py" masif_site "$1"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/04-masif_precompute.py" masif_ppi_search "$1"
fi


