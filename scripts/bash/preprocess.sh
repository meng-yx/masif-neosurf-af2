#!/bin/bash

PDB_PATH=$1
TARGET=$2
LIGAND=$3
LIGAND_PATH=$4
OUTPUT_DIR=$5

CONFIG_FILE="${CONFIG_FILE:-scripts/config.sh}"
if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "Error: run from masif-neosurf repo root (${CONFIG_FILE} not found; run from repo root, set CONFIG_FILE to override)" >&2
    exit 1
fi

source "${CONFIG_FILE}"

PREPROCESS_PY_ARGS=(--infer_reduce_het_dict)
if [[ "${PREPROCESS_OVERWRITE:-0}" == "1" ]]; then
    PREPROCESS_PY_ARGS+=(--overwrite)
fi

apptainer exec --writable-tmpfs -B "${BIND_MOUNTS}" --env "REDUCE_HET_DICT=${REDUCE_HET_DICT}" "${IMAGE}" \
    python -W ignore preprocess_pdb.py \
        "${PDB_PATH}" \
        "${TARGET}" \
        -l "${LIGAND}" \
        -s "${LIGAND_PATH}" \
        -o "${OUTPUT_DIR}" \
        "${PREPROCESS_PY_ARGS[@]}"
