#!/usr/bin/env bash
# Usage:
#   ./train_nn.sh nn_models.v1.custom_params
#
# Runs training inside the conda env defined by environment.yml in this
# directory (default name: masif-coherence). Create it first:
#   conda env create -f environment.yml
#
# Override env name if needed:
#   CONDA_ENV_NAME=my-env ./train_nn.sh

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
    echo "WARNING: conda not found; using python3 from PATH (install deps via pip -r .../requirements.txt)" >&2
    python3 "$@"
  fi
}

_run_python "$masif_seed_search_source/masif_coherence/train.py" "${1:-nn_models.v1.custom_params}"
