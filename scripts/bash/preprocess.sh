#!/bin/bash

PDB_PATH=$1
TARGET=$2
LIGAND=$3
LIGAND_PATH=$4
OUTPUT_DIR=$5

if [[ ! -f scripts/config.sh ]]; then
    echo "Error: run from masif-neosurf repo root (scripts/config.sh not found in $(pwd))" >&2
    exit 1
fi

source scripts/config.sh

srun apptainer exec --writable-tmpfs -B "${BIND_MOUNTS}" --env "REDUCE_HET_DICT=${REDUCE_HET_DICT}" "${IMAGE}" \
    python -W ignore preprocess_pdb.py \
        "${PDB_PATH}" \
        "${TARGET}" \
        -l "${LIGAND}" \
        -s "${LIGAND_PATH}" \
        -o "${OUTPUT_DIR}"
