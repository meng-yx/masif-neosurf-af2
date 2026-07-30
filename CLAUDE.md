# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

This is a fork of **MaSIF-neosurf** (LPDI-EPFL), the surface-based deep-learning tool from Marchand et al., *Nature* 2025 ("Targeting protein-ligand neosurfaces with a generalizable deep learning tool"). The upstream tool preprocesses protein-ligand complexes into learned molecular surface fingerprints and searches for complementary binding "seeds".

On top of that upstream base, this fork adds a **large-scale "Neosurf-on-Neosurf" screen**: an all-vs-all search where many ligand-bound "target" complexes are queried against a large library of ligand-bound "seed" complexes (e.g. the human reference proteome). Every predicted interaction pairs a different target *and* a different seed, each carrying its own ligand. This new pipeline lives in `scripts/`, `notebooks/`, and the config files; the upstream MaSIF code (`masif/`, `masif_seed_search/`, `masif_search.py`, `preprocess_pdb.py`, `rosetta_scripts/`, `computational_benchmark/`) is largely inherited.

When working here, assume the large-scale screen (`scripts/` + `notebooks/neosurf_neosurf.ipynb`) is the active work unless the task clearly concerns upstream single-target search.

## Two execution environments — do not confuse them

The upstream MaSIF code (preprocessing, surface computation, search) runs **only inside a container** because it depends on MSMS, APBS, PyMesh, reduce, and TensorFlow 1.x. This fork runs that container via **Apptainer/Singularity** (`masif-neosurf_v0.1.sif`), not Docker:

```bash
apptainer exec -B /work:/work,/scratch:/scratch masif-neosurf_v0.1.sif python -W ignore masif_search.py ...
```

The orchestration layer (the `scripts/python/*.py` manifest drivers, the notebook, the Dash app) runs **outside** the container in a conda env named `MaSIF` (`conda activate MaSIF`), and shells into the container for the heavy MaSIF steps. The SLURM scripts do both: `conda activate MaSIF` for the Python driver, which then calls `apptainer exec` per row.

The upstream README documents a Docker workflow (`docker build . -t masif-neosurf`); that is the vanilla single-target path and still valid, but the large-scale scripts assume the Apptainer `.sif` image named in `scripts/config.sh`.

## The large-scale pipeline (this fork's core)

Orchestrated end-to-end from `notebooks/neosurf_neosurf.ipynb`. Four numbered stages, each fanned out as a SLURM array job over round-robin CSV subsets (`input_{k}.csv` -> `output_{k}.csv`, one file per array task). Intermediate manifests accumulate under `data/processing/{1_prep_input,2_masif_preprocess,3_masif_search}/`.

1. **Prepare input** — `scripts/python/prepare_input.py` (via `scripts/slurm/prepare_input_array.sh`). For each `(pdb_id, protein_chain, ligand_chain, ligand_code)` row: downloads structure + ligand SDF from RCSB, trims to target chain + one ligand residue, runs **EvoEF2 RepairStructure**, and merges the ligand back onto the repaired protein. Emits `.pdb`/`.sdf` into `data/input/` plus manifest columns. **Read this script's module docstring before touching identifier logic** — it defines the identifier scheme (`target = {pdb_id}-{ligand_code}_{pdb_chain_suffix}` is the canonical `seed_id`/`ppi_pair_id`; author chain IDs vs PDB-normalized chain IDs; long ligand codes stored as `ligand_resname`). Getting this wrong silently breaks the join between search results and metadata.

2. **MaSIF preprocess** — `scripts/python/preprocess_manifest.py` (via `preprocess_array.sh`), which shells into `scripts/bash/preprocess.sh` -> containerized `preprocess_pdb.py`. Produces surface descriptors under `data/preprocess/`. Idempotent: rows are **skipped** when `descriptors/sc05/all_feat/{target}/p1_desc_{straight,flipped}.npy` already exist. Each row is appended to the output CSV immediately (crash-safe resume).

3. **MaSIF search** — `scripts/bash/search.sh` (via `search_array.sh`) -> containerized `masif_search.py`. This is the all-vs-all step. `search_array.sh` has two modes: single-target (target id as arg 1) or **multi-target** (a `.txt` manifest of targets as arg 1, looping every target per array task). Output lands under `data/masif_search/{target}/`, with per-anchor-site `site_*/` dirs and final `clustered_matches/*.csv`.

4. **Analysis** — the notebook gathers all `data/masif_search/*/clustered_matches/*.csv` into `df_results_all.csv`, joins metadata from the preprocess manifest (`df_preprocess_ok.csv`), and deduplicates to `df_results_dedup.csv` (one row per `target`×`matched_protein`×`cluster_id`, keeping highest score).

### Neosurf-specific search flags (this fork's key additions to `masif_search.py`)

- `--target-auto-neosurf` / `--seed-auto-neosurf`: automatically anchor the search site on the **largest HETATM (ligand) residue** instead of a manually specified `--target_chain/--target_resid/--target_atom_id`. This is what makes the screen "neosurf-on-neosurf". Controlled by `TARGET_AUTO_NEOSURF` / `SEED_AUTO_NEOSURF` in config.
- `--seed_subset <file>`: restrict a search task to the seed ids listed in one subset file (how the array fan-out shards the seed library).
- `--resume`: per-site completion tracking so re-submitted array jobs skip finished work (uses `search_output.is_site_completed` / `mark_site_completed`, not a pile of empty CSVs).

## Configuration

All pipeline parameters are centralized in **`scripts/config.sh`**, sourced by every SLURM/bash script from the repo root (scripts hard-fail if `scripts/config.sh` is not found in `$(pwd)`, i.e. **always launch from repo root**). Key knobs: MaSIF-search cutoffs (`CUTOFF`, `NN_SCORE_CUTOFF`, `RANSAC_ITER`, `TARGET_SAMPLING_RATIO`), the auto-neosurf toggles, `RESUME`, the Apptainer `IMAGE`, `BIND_MOUNTS`, and the `DATABASE_ROOT`/`TARGET_PREPROCESS_ROOT` (both `data/preprocess`). `scripts/configs/config_*.sh` are named presets (currently identical copies) for different screens.

## Commands

```bash
# One-time: build EvoEF2 (needed by stage 1 prepare_input)
bash scripts/slurm/install_dependencies.sh

# Launch a stage (always from repo root). N = number of round-robin subsets.
sbatch --array=1-N scripts/slurm/prepare_input_array.sh \
    data/processing/1_prep_input/input_subsets data/input \
    data/processing/1_prep_input/output_subsets EvoEF2/EvoEF2
sbatch --array=1-N scripts/slurm/preprocess_array.sh \
    data/processing/2_masif_preprocess/input_subsets \
    data/processing/2_masif_preprocess/output_subsets
sbatch --array=1-N scripts/slurm/search_array.sh \
    data/masif_search/query_targets.txt data/masif_search data/masif_search/subset

# Run one stage locally without SLURM (array id defaults to 1)
python scripts/python/preprocess_manifest.py --input_csv ... --out_csv ... --preprocess_root data/preprocess --repo_root "$(pwd)"

# Vanilla upstream single-complex search (inside container)
apptainer exec -B /work:/work,/scratch:/scratch masif-neosurf_v0.1.sif \
    python -W ignore masif_search.py --target 8VLB_A --target_dir example/processed ...

# Interactive Jupyter on a compute node (prints a port to tunnel)
sbatch scripts/slurm/jupyter.sh

# Restore gitignored data (data/, .sif, EvoEF2/) from the /work backup into this scratch clone
sbatch scripts/slurm/restore_from_work.sh            # add --dry-run to preview
```

There is no test suite, linter, or build step for the Python orchestration scripts — validation is done by running a small subset (the notebook uses `N_SEED` / a single `pdb_id` filter to smoke-test one complex before scaling up).

## Results explorer app

`scripts/python/neosurf_neosurf_explorer_app.py` is a **Dash** app (run `python scripts/python/neosurf_neosurf_explorer_app.py`, serves on the default Dash port with `debug=True`). It reads `df_preprocess_ok.csv` (per-structure metadata keyed by `target`), `df_results_dedup.csv` (one row per cluster — drives the tree + scatter), and `df_results_all.csv` (every pose, loaded once and indexed by cluster key, populated lazily on cluster expand). It renders a lazy collapsible tree (target-protein → target → matched-protein → seed → cluster → pose), a Plotly scatter, and a py3Dmol viewer showing target + transform-superposed seed with both ligands as sticks. `pdb_path` values are rewritten at load time to `data/preprocess/data_preparation/01-benchmark_pdbs/`. The module docstring is the authoritative spec for the data model — read it before changing data loading.

`scripts/python/results_to_pdb.py` reproduces the app's superposition offline: it applies each result row's flattened 4×4 transform to the matched protein and writes assembled target+seed complexes as `.pdb` files (BioPython, mirroring the PyMOL command sequence in `notebooks/scratch.ipynb`).

## Data layout

Everything under `data/` (and `EvoEF2/`, `*.sif`, `logs/`) is **gitignored** and lives on scratch, restorable from a `/work` backup via `restore_from_work.sh`. Notable roots: `data/input/` (prepared PDBs+SDFs), `data/preprocess/` (MaSIF surfaces + descriptors, the search database), `data/processing/{1,2,3}_*/` (per-stage manifests + round-robin subsets), `data/masif_search/{target}/clustered_matches/` (final hits).

## Gotchas

- **Always run scripts from the repo root** — they source `scripts/config.sh` by relative path and exit otherwise.
- The join between search results and metadata is keyed on the `target` id, which equals the `seed_id`/`ppi_pair_id` computed by `prepare_input.py`. Never re-derive `pdb_id` or chain from the `target` string by splitting on `-`/`_` — ligand codes and multi-letter author chains make that ambiguous. Use the manifest lookup.
- Ligand preprocessing is the fragile part (documented in the README "Known limitations"): `reduce` needs an up-to-date `REDUCE_HET_DICT` (set in config to a `/work` path), some ligands need a hand-edited mol2 passed via `-m`, and ring-heavy ligands can break pdb2pqr (patch at `masif/data/masif_neosurf/pdb2pqr_mol2_patch.py`).
- Stages 2 and 3 are resumable/idempotent by design (descriptor-existence check; `--resume` site markers). Re-submitting an array job over the same subsets continues rather than redoing work.
