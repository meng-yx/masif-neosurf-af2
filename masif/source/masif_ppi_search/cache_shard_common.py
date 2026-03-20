import hashlib
import importlib
import json
import os
import sys

import numpy as np

from default_config.masif_opts import masif_opts


FINAL_CACHE_FILES = [
    "pos_names.npy",
    "neg_names.npy",
    "binder_rho_wrt_center.npy",
    "binder_theta_wrt_center.npy",
    "binder_input_feat.npy",
    "binder_mask.npy",
    "pos_training_idx.npy",
    "pos_val_idx.npy",
    "pos_test_idx.npy",
    "pos_rho_wrt_center.npy",
    "pos_theta_wrt_center.npy",
    "pos_input_feat.npy",
    "pos_mask.npy",
    "neg_training_idx.npy",
    "neg_val_idx.npy",
    "neg_test_idx.npy",
    "neg_rho_wrt_center.npy",
    "neg_theta_wrt_center.npy",
    "neg_input_feat.npy",
    "neg_mask.npy",
]


def parse_cli(custom_params_required=False):
    cli = {
        "custom_params_module": None,
        "extra_args": [],
    }
    if len(sys.argv) > 1 and not sys.argv[1].startswith("-"):
        cli["custom_params_module"] = sys.argv[1]
        cli["extra_args"] = sys.argv[2:]
    else:
        cli["extra_args"] = sys.argv[1:]
    if custom_params_required and cli["custom_params_module"] is None:
        raise ValueError("A custom params module is required as first positional arg.")
    return cli


def load_params(custom_params_module=None):
    params = dict(masif_opts["ppi_search"])
    if custom_params_module is not None:
        custom_params = importlib.import_module(custom_params_module, package=None)
        custom_params = custom_params.custom_params
        for key in custom_params:
            print("Setting {} to {} ".format(key, custom_params[key]))
            params[key] = custom_params[key]
    if "pids" not in params:
        params["pids"] = ["p1", "p2"]
    return params


def normalize_mix(w_cross, w_within, w_hard):
    mix = np.asarray([w_cross, w_within, w_hard], dtype=float)
    mix = np.maximum(mix, 0.0)
    if mix.sum() <= 0.0:
        mix = np.asarray([0.0, 1.0, 0.0], dtype=float)
    mix /= mix.sum()
    return mix


def split_counts(total, mix):
    raw = total * mix
    counts = np.floor(raw).astype(int)
    remainder = int(total - counts.sum())
    if remainder > 0:
        order = np.argsort(-(raw - counts))
        counts[order[:remainder]] += 1
    return counts


def sample_from_index_array(index_array, n_samples):
    if n_samples <= 0 or len(index_array) == 0:
        return np.asarray([], dtype=int)
    replace = len(index_array) < n_samples
    return np.random.choice(index_array, size=n_samples, replace=replace)


def concat_or_empty(chunks, ndim):
    if len(chunks) == 0:
        if ndim == 2:
            return np.empty((0, 0))
        if ndim == 3:
            return np.empty((0, 0, 0))
        return np.empty((0,))
    return np.concatenate(chunks, axis=0)


def remap_indices(indices, keep_mask):
    indices = np.asarray(indices, dtype=int)
    if len(indices) == 0:
        return indices
    old_to_new = -1 * np.ones(len(keep_mask), dtype=int)
    old_to_new[np.where(keep_mask)[0]] = np.arange(np.sum(keep_mask))
    valid_old = keep_mask[indices]
    remapped = old_to_new[indices[valid_old]]
    assert not (remapped == -1).any()
    return remapped


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def to_native(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, dict):
        return {k: to_native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_native(v) for v in value]
    return value


def params_fingerprint(payload):
    normalized = to_native(payload)
    text = json.dumps(normalized, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def save_json(path, payload):
    with open(path, "w") as handle:
        json.dump(to_native(payload), handle, indent=2, sort_keys=True)


def read_json(path):
    with open(path, "r") as handle:
        return json.load(handle)


def write_success_marker(directory):
    with open(os.path.join(directory, "_SUCCESS"), "w") as handle:
        handle.write("ok\n")
