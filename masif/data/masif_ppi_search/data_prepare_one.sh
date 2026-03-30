masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_root=$masif_neosurf_root/masif
masif_source=$masif_root/source/
docker_image=$masif_neosurf_root/masif-neosurf_v0.1.sif
#echo "docker image: $docker_image"
export PYTHONPATH=$PYTHONPATH:$masif_source
if [ -z "$2" ]; then
    echo "Usage: $0 PDBID_CHAIN1_CHAIN2 /path/to/pinder_pdb_dir"
    exit 1
fi
MODEL=$(echo $1| cut -d"_" -f1)
CHAIN1=$(echo $1| cut -d"_" -f2)
CHAIN2=$(echo $1| cut -d"_" -f3)
PDB_ID=$(echo $MODEL | cut -d"-" -f1)
PDB_CHAIN_L=$(echo $MODEL| cut -d"-" -f2)
PDB_CHAIN_R=$(echo $MODEL| cut -d"-" -f3)
MODEL_TYPE=$(echo $MODEL| cut -d"-" -f4)
PINDER_PDB_DIR=$2
PINDER_PDB_DIR_REAL=$(readlink -f "$PINDER_PDB_DIR")

if [ ! -d "$PINDER_PDB_DIR_REAL" ]; then
    echo "Error: pinder_pdb_dir does not exist or is not a directory: $PINDER_PDB_DIR"
    exit 1
fi

# Bind repo root and the resolved Pinder directory (important when pinderMaSIF is a symlink).
SINGULARITY_BIND="$masif_neosurf_root:$masif_neosurf_root,$PINDER_PDB_DIR_REAL:$PINDER_PDB_DIR_REAL"

run_prepare_pipeline() {
    local COMPLEX_ID="$1"
    local COMPLEX_MODEL
    local COMPLEX_CHAIN1
    local COMPLEX_CHAIN2

    COMPLEX_MODEL=$(echo "$COMPLEX_ID" | cut -d"_" -f1)
    COMPLEX_CHAIN1=$(echo "$COMPLEX_ID" | cut -d"_" -f2)
    COMPLEX_CHAIN2=$(echo "$COMPLEX_ID" | cut -d"_" -f3)

    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/00d-pdb_from_pinder.py" "$COMPLEX_ID" --pinder-pdb-dir "$PINDER_PDB_DIR_REAL"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/01-pdb_extract_and_triangulate.py" "${COMPLEX_MODEL}_${COMPLEX_CHAIN1}"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/01-pdb_extract_and_triangulate.py" "${COMPLEX_MODEL}_${COMPLEX_CHAIN2}"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/04-masif_precompute.py" masif_site "$COMPLEX_ID"
    singularity exec --bind "$SINGULARITY_BIND" "$docker_image" python "$masif_source/data_preparation/04-masif_precompute.py" masif_ppi_search "$COMPLEX_ID"
}

has_precompute_outputs() {
    local COMPLEX_ID="$1"
    local COMPLEX_MODEL
    local COMPLEX_CHAIN1
    local COMPLEX_CHAIN2
    local OUT1
    local OUT2

    COMPLEX_MODEL=$(echo "$COMPLEX_ID" | cut -d"_" -f1)
    COMPLEX_CHAIN1=$(echo "$COMPLEX_ID" | cut -d"_" -f2)
    COMPLEX_CHAIN2=$(echo "$COMPLEX_ID" | cut -d"_" -f3)

    OUT1="data_preparation/04b-precomputation_12A/precomputation/${COMPLEX_MODEL}_${COMPLEX_CHAIN1}_${COMPLEX_CHAIN2}/p1_sc_labels.npy"
    OUT2="data_preparation/04b-precomputation_12A/precomputation/${COMPLEX_MODEL}_${COMPLEX_CHAIN1}_${COMPLEX_CHAIN2}/p2_sc_labels.npy"
    if [ -f "$OUT1" ] && [ -f "$OUT2" ]; then
        echo "$OUT1 and $OUT2 already exist"
        return 0
    fi
    return 1
}

#echo "Using Pinder PDB directory: $PINDER_PDB_DIR_REAL"


# Pinder-specific logic for handling aligned AF2 models (AP/PA)
# ---- New: check if the complex_id is a PP model -----
if [ "$MODEL_TYPE" == "PP" ]; then
    if has_precompute_outputs "$1"; then
        echo "Skipped PP data preparation: ${1}"
        exit 0
    fi
    echo "Preparing PP model: ${1}"
    run_prepare_pipeline "$1"
    exit 0
else
    
    PP_MODEL="${PDB_ID}-${PDB_CHAIN_L}-${PDB_CHAIN_R}-PP_L_R"
   
    if has_precompute_outputs "$PP_MODEL"; then
        echo "Skipped equivalent PP model data preparation: ${PP_MODEL}"
    else
        echo "Preparing equivalent PP model: ${PP_MODEL}"
        run_prepare_pipeline "${PP_MODEL}"
    fi
    if has_precompute_outputs "$1"; then
        echo "Skipped target model data preparation: ${1}"
    else
        echo "Preparing target model: ${1}"
        run_prepare_pipeline "$1"
    fi
fi


