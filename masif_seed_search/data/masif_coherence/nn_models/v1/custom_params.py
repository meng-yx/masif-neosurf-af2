"""v1 baseline experiment overrides for MaSIF-Coherence.

The shell runner (`train_nn.sh`) passes this module's path to `train.py`, which
imports it and overlays `custom_params` onto the defaults in
`masif_seed_search/source/masif_coherence/config.py`.
"""

custom_params = {
    "model_dir": "nn_models/v1/model_data/",
    # The defaults already encode the v1 baseline; this file exists so downstream
    # experiments can be started by copying this directory and editing only the
    # fields that differ.
}
