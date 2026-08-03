#!/usr/bin/env bash
# Run masif_search for one SLURM array task using a seed subset file.
#
# Usage:
#   bash scripts/bash/search.sh <query_target> <out_dir> <seed_subset_file>
#
# When launched from search_array.sh, seed_subset_file is
#   ${seed_subset_dir}/${SLURM_ARRAY_TASK_ID}

set -euo pipefail

QUERY_TARGET=$1
OUT_DIR=$2
SEED_SUBSET=$3

if [[ -z "${QUERY_TARGET}" || -z "${OUT_DIR}" || -z "${SEED_SUBSET}" ]]; then
    echo "Usage: bash scripts/bash/search.sh <query_target> <out_dir> <seed_subset_file>" >&2
    exit 1
fi

CONFIG_FILE="${CONFIG_FILE:-scripts/config.sh}"
if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "Error: run from masif-neosurf repo root (${CONFIG_FILE} not found; run from repo root, set CONFIG_FILE to override)" >&2
    exit 1
fi

if [[ ! -f "${SEED_SUBSET}" ]]; then
    echo "Error: seed subset file not found: ${SEED_SUBSET}" >&2
    exit 1
fi

source "${CONFIG_FILE}"

AUTO_FLAGS=()
if [[ "${TARGET_AUTO_NEOSURF:-0}" == "1" ]]; then
    AUTO_FLAGS+=(--target-auto-neosurf)
fi
if [[ "${SEED_AUTO_NEOSURF:-0}" == "1" ]]; then
    AUTO_FLAGS+=(--seed-auto-neosurf)
    AUTO_FLAGS+=(--seed_site_cutoff "${CUTOFF}")
fi
if [[ "${RESUME:-0}" == "1" ]]; then
    AUTO_FLAGS+=(--resume)
fi

echo "Search ${QUERY_TARGET} seeds from ${SEED_SUBSET}"

apptainer exec -B "${BIND_MOUNTS}" "${IMAGE}" \
    python -W ignore masif_search.py \
        --target "${QUERY_TARGET}" \
        --target_dir "${TARGET_PREPROCESS_ROOT}" \
        --out_dir "${OUT_DIR}" \
        --target_site_cutoff "${CUTOFF}" \
        --target_site_sample_ratio "${TARGET_SAMPLING_RATIO}" \
        --seed_dir "${DATABASE_ROOT}" \
        --seed_subset "${SEED_SUBSET}" \
        --seed_iface_cutoff "${IFACE_CUTOFF}" \
        --seed_desc_dist_cutoff "${DESC_DIST_CUTOFF}" \
        --seed_nn_score_cutoff "${NN_SCORE_CUTOFF}" \
        --ransac_iter "${RANSAC_ITER}" \
        --n_retry_alignment "${N_RETRY_ALIGNMENT}" \
        "${AUTO_FLAGS[@]}"
