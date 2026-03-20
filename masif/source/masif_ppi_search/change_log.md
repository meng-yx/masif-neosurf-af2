# MaSIF PPI Cache Refactor Change Log

This document summarizes the refactor from the original single-process
`masif_ppi_search_cache_training_data.py` workflow to the current 3-step
pipeline in `masif/source/masif_ppi_search`.

## Why this refactor was done

The original cache generation ran in one process and performed all steps in a
single pass:

1. discover valid PPI pairs
2. assign train/val/test split
3. build positive and negative candidates
4. sample negatives
5. accumulate all tensors in memory
6. concatenate and write final cache arrays

At dataset scale, this caused very large memory peaks and frequent OOM kills,
especially during global pooling and end-of-job concatenation.

## Old architecture (single process)

- Script: `masif_ppi_search_cache_training_data.py` (monolithic mode)
- Characteristics:
  - one job performs both metadata building and tensor materialization
  - global pools and sampled tensors are retained in process memory
  - final `np.concatenate(...)` over large lists creates additional peak RAM
- Output:
  - final training-compatible cache arrays directly in `params["cache_dir"]`

## New architecture (3-step workflow)

The workflow is now split into three scripts with separate execution stages:

1. **Catalog build**: `build_cache_catalog.py`
2. **Shard generation**: `masif_ppi_search_cache_training_data.py`
3. **Shard merge**: `merge_cache_shards.py`

### Stage 1: Catalog build (`build_cache_catalog.py`)

Responsibilities:

- scans precomputation pairs and validates usable examples
- computes record-level metadata:
  - `split` (`train` / `val` / `test`)
  - positive index pairs (`k1`, `k2`)
  - within-complex negative candidates (`k_neg2`, `within_dneg`)
- writes compact catalog artifacts under `<cache_dir>/catalog`

Refactor details:

- pool files were compacted to numeric mappings only:
  - `*_source_record_idx.npy` (int32)
  - `*_local_within_idx.npy` (int32)
- object-heavy arrays (`*_ppi.npy`, `*_pdb.npy`, `*_names.npy`) were removed
- catalog manifest now fingerprints key parameters and acts as contract input

### Stage 2: Shard generation (`masif_ppi_search_cache_training_data.py`)

Responsibilities:

- loads catalog and processes only records assigned to one shard:
  - assignment rule: `record_idx % num_subsets == subset_id`
- preserves original negative sampling semantics:
  - `split_counts(neg_ratio, mix)`
  - cross / within / hard buckets
  - fallback order: `cross_valid -> within -> cross_all`
- writes subset artifacts in shard-local directory:
  - `<subset_root>/subset_<id>/parts/...`
  - `<subset_root>/subset_<id>/manifest.json`
  - `<subset_root>/subset_<id>/_SUCCESS`

Memory-oriented refactor details:

- numeric catalog arrays are loaded with `mmap_mode='r'`
- metadata such as cross names are derived lazily from record indices
- unbounded pair cache replaced by bounded LRU cache (`PairArrayCache`)
- output writing changed from one-shot concatenation to chunked part writes
  (`SubsetPartWriter`) with periodic flushes
- NaN filtering and index remapping happen per part before persistence

### Stage 3: Merge (`merge_cache_shards.py`)

Responsibilities:

- validates every expected subset and `_SUCCESS` marker
- verifies catalog fingerprint consistency across subset manifests
- loads and concatenates subset outputs in deterministic subset order
- rebases index arrays with cumulative positive/negative offsets
- writes final cache arrays in exactly the format expected by training

Compatibility:

- merge supports both:
  - new chunked subset format (`format_version >= 2`)
  - legacy flat subset files (`format_version == 1`)

## Data contract and downstream compatibility

Final outputs remain compatible with `masif_ppi_search_train.py`, including:

- binder arrays: `binder_*`
- positive arrays: `pos_*`
- negative arrays: `neg_*`
- split indices: `pos_*_idx.npy`, `neg_*_idx.npy`
- names: `pos_names.npy`, `neg_names.npy`

This preserves downstream training behavior without requiring training-script
changes.

## Slurm and wrappers

Separate entrypoints were added for each stage:

- catalog:
  - `masif/data/masif_ppi_search/build_cache_catalog.sh`
  - `masif/data/masif_ppi_search/build_cache_catalog.slurm`
- shard array:
  - `masif/data/masif_ppi_search/cache_nn_shards.sh`
  - `masif/data/masif_ppi_search/cache_nn_shards.slurm`
- merge:
  - `masif/data/masif_ppi_search/merge_cache_shards.sh`
  - `masif/data/masif_ppi_search/merge_cache_shards.slurm`

Runtime knobs introduced for shard workers:

- `NUM_SUBSETS`
- `PAIR_CACHE_MAX_RECORDS`
- `SHARD_FLUSH_EVERY_POS`

These knobs allow tuning memory footprint vs. I/O throughput.

## Behavioral guarantees kept intact

The refactor intentionally preserves:

- split logic from training/testing lists and `range_val_samples`
- positive contact extraction logic
- negative sampling mix and fallback semantics
- split disjointness assertions
- final cache file naming and structure

## Net impact

- significantly lower per-process memory pressure in shard jobs
- better horizontal scaling through Slurm arrays
- improved fault isolation (retry individual subsets)
- deterministic merge point with contract checks before final cache publication
