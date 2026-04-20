# MaSIF-Coherence

A protein-level interaction scorer built on top of `masif_ppi_search` descriptors. Instead of scoring a rigid alignment between two patches (MaSIF-search's second stage), MaSIF-Coherence reasons relationally over a set of candidate vertex-pair matches and learns whether those matches are mutually geometrically consistent.

See:

- `masif_coherence_design.md` — full design rationale and SE(3)-invariant architecture.
- `minimal_v1_plan.md` — minimal v1 baseline this codebase implements.

## Layout

- Code: `masif_seed_search/source/masif_coherence/` (PyTorch).
- Per-experiment overrides: `nn_models/<exp>/custom_params.py`.
- Per-experiment checkpoints: `nn_models/<exp>/model_data/`.
- Pair lists: `lists/training.txt`, `lists/testing.txt` (re-used from `masif/data/masif_ppi_search/lists/` in v1).

## Environment

Create the conda environment from this directory (same layout as other `train_nn.sh` wrappers, but PyTorch instead of the TF1 Singularity image):

```
conda env create -f masif_seed_search/data/masif_coherence/environment.yml
```

`train_nn.sh` runs `python` via `conda run -n masif-coherence` so you do not need to `conda activate` first. Override the env name with `CONDA_ENV_NAME=... ./train_nn.sh` if needed.

## Usage

Train v1 on the default experiment:

```
cd masif_seed_search/data/masif_coherence
./train_nn.sh nn_models.v1.custom_params
```

Sanity-check SE(3) and permutation invariance before any training:

```
conda run -n masif-coherence --no-capture-output python masif_seed_search/source/masif_coherence/check_invariance.py
```

## Requirements

Preferred: conda (see above). Alternative pip-only install:

```
pip install -r masif_seed_search/source/masif_coherence/requirements.txt
```

This module is self-contained and does not modify any existing MaSIF (TF1) environment.
