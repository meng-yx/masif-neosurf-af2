#!/usr/bin/env bash
# Usage:
#   ./build_cache.sh nn_models.v1.custom_params

set -euo pipefail

masif_neosurf_root=$(git rev-parse --show-toplevel)
masif_seed_search_source=$masif_neosurf_root/masif_seed_search/source
coherence_data_root=$masif_neosurf_root/masif_seed_search/data/masif_coherence
CONDA_ENV_NAME="${CONDA_ENV_NAME:-masif-coherence}"

export PYTHONPATH=$masif_seed_search_source:$coherence_data_root:${PYTHONPATH:-}

cd "$coherence_data_root"

_run_python() {
  if command -v conda >/dev/null 2>&1; then
    conda run -n "$CONDA_ENV_NAME" --no-capture-output python "$@"
  else
    echo "WARNING: conda not found; using python3 from PATH." >&2
    python3 "$@"
  fi
}

_run_python "$masif_seed_search_source/masif_coherence/build_cache.py" "${1:-nn_models.v1.custom_params}"
