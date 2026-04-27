"""v1 baseline experiment overrides for MaSIF-Coherence.

The shell runner (`train_nn.sh`) passes this module's path to `train.py`, which
imports it and overlays `custom_params` onto the defaults in
`masif_seed_search/source/masif_coherence/config.py`.
"""

custom_params = {
    # ------------------------------------------------------------------
    # Data paths (relative to the directory train.py is run from, matching
    # the masif_opts["ppi_search"] convention).
    # ------------------------------------------------------------------
    "descriptors_dir": "descriptors/sc05/all_feat/",
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
