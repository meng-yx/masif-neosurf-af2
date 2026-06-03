#!/usr/bin/env bash
# Run masif_search for one row of the search manifest.
#
# Usage:
#   bash scripts/bash/search.sh \
#     <query_target> <target_chain> <target_resid> \
#     <seed_target> <seed_chain> <seed_resid> <out_dir>

set -euo pipefail

QUERY_TARGET=$1
TARGET_CHAIN=$2
TARGET_RESID=$3
SEED_TARGET=$4
SEED_CHAIN=$5
SEED_RESID=$6
OUT_DIR=$7

if [[ ! -f scripts/config.sh ]]; then
    echo "Error: run from masif-neosurf repo root (scripts/config.sh not found in $(pwd))" >&2
    exit 1
fi

source scripts/config.sh

mkdir -p "${OUT_DIR}"
SEED_SUBSET=$(mktemp)
trap 'rm -f "${SEED_SUBSET}"' EXIT
printf '%s\n' "${SEED_TARGET}" > "${SEED_SUBSET}"

srun apptainer exec -B "${BIND_MOUNTS}" "${IMAGE}" \
    python -W ignore masif_search.py \
        --target "${QUERY_TARGET}" \
        --target_dir "${DATABASE_ROOT}" \
        --out_dir "${OUT_DIR}" \
        --target_chain "${TARGET_CHAIN}" \
        --target_resid "${TARGET_RESID}" \
        --target_site_cutoff "${CUTOFF}" \
        --target_site_sample_ratio "${TARGET_SAMPLING_RATIO}" \
        --seed_dir "${DATABASE_ROOT}" \
        --seed_subset "${SEED_SUBSET}" \
        --seed_chain "${SEED_CHAIN}" \
        --seed_resid "${SEED_RESID}" \
        --seed_site_cutoff "${CUTOFF}" \
        --seed_iface_cutoff "${IFACE_CUTOFF}" \
        --seed_desc_dist_cutoff "${DESC_DIST_CUTOFF}" \
        --seed_nn_score_cutoff "${NN_SCORE_CUTOFF}" \
        --ransac_iter "${RANSAC_ITER}" \
        --n_retry_alignment "${N_RETRY_ALIGNMENT}"
