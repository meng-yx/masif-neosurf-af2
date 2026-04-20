"""Default hyperparameters for MaSIF-Coherence v1.

These are the baseline values committed to by `minimal_v1_plan.md`. Per-experiment
overrides live next to each checkpoint in
`masif_seed_search/data/masif_coherence/nn_models/<exp>/custom_params.py`, mirroring
the `custom_params.py` convention used by `masif/data/masif_ppi_search/nn_models/sc05/`.
"""

from copy import deepcopy


coherence_opts = {
    # ------------------------------------------------------------------
    # Data paths (relative to the directory train.py is run from, matching
    # the masif_opts["ppi_search"] convention).
    # ------------------------------------------------------------------
    "masif_precomputation_dir": "data_preparation/04b-precomputation_12A/precomputation/",
    "ply_chain_dir": "data_preparation/01-benchmark_surfaces/",
    "training_list": "lists/training.txt",
    "testing_list": "lists/testing.txt",
    "model_dir": "nn_models/v1/model_data/",
    # ------------------------------------------------------------------
    # Model hyperparameters.
    # ------------------------------------------------------------------
    "descriptor_dim": 80,
    "top_k": 256,
    "hidden_dim": 64,
    "num_heads": 4,
    "num_layers": 2,
    "ffn_mult": 2,
    "rbf_num_bins": 16,
    "rbf_min_ang": 0.0,
    "rbf_max_ang": 40.0,
    # ------------------------------------------------------------------
    # Training hyperparameters.
    # ------------------------------------------------------------------
    "seed": 42,
    "learning_rate": 1e-4,
    "weight_decay": 0.0,
    "grad_accum_steps": 16,
    "num_epochs": 50,
    "val_fraction": 0.1,
    "log_every": 50,
    "save_every_epoch": True,
    # One random non-cognate pair per positive, per the v1 plan.
    "neg_per_pos": 1,
    # Device: "cuda" if available, else "cpu".
    "device": "auto",
}


def get_default_params():
    """Return a fresh deep-copy of the default v1 parameters."""
    return deepcopy(coherence_opts)


def merge_custom_params(base, override):
    """Recursively overlay `override` onto `base` (mirrors the helper in
    `masif_ppi_search_train.py`). Returns the merged dict for convenience.
    """
    for key, value in override.items():
        if key in base and isinstance(base[key], dict) and isinstance(value, dict):
            merge_custom_params(base[key], value)
        else:
            if key not in base:
                print(f"Adding new key {key}={value}.")
            else:
                print(f"Setting {key} to {value}.")
            base[key] = value
    return base
