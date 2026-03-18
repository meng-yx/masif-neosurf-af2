# Header variables and parameters.
import importlib
import os
import sys

import numpy as np
import pymesh
from scipy.spatial import cKDTree

from default_config.masif_opts import masif_opts

"""
masif_ppi_search_cache_training_data.py: Function to cache all the training data for MaSIF-search. 
                This function extract all the positive pairs and a random number of negative surfaces.
                In the future, the number of negative surfaces should be increased.
Pablo Gainza - LPDI STI EPFL 2019
Released under an Apache License 2.0

Modification: Yanxiang Meng - 18.3.2026
Support configurable negative sampling with mixed sources:
within-complex, cross-complex, and hard negatives (near-interface).
"""
def normalize_mix(w_cross, w_within, w_hard):
    mix = np.asarray([w_cross, w_within, w_hard], dtype=float)
    mix = np.maximum(mix, 0.0)
    if mix.sum() <= 0.0:
        # Backward-compatible default: all negatives from within-complex.
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


def sample_indices(pool_size, n_samples):
    if n_samples <= 0 or pool_size <= 0:
        return np.asarray([], dtype=int)
    replace = pool_size < n_samples
    return np.random.choice(pool_size, size=n_samples, replace=replace)


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


params = masif_opts["ppi_search"]
if len(sys.argv) > 1:
    custom_params_file = sys.argv[1]
    custom_params = importlib.import_module(custom_params_file, package=None)
    custom_params = custom_params.custom_params
    for key in custom_params:
        print("Setting {} to {} ".format(key, custom_params[key]))
        params[key] = custom_params[key]

if "pids" not in params:
    params["pids"] = ["p1", "p2"]

neg_ratio = int(params.get("neg_ratio", 1))
neg_ratio = max(1, neg_ratio)
mix = normalize_mix(
    params.get("neg_mix_cross_complex", 0.0),
    params.get("neg_mix_within_complex", 1.0),
    params.get("neg_mix_hard", 0.0),
)
enforce_diff_pdb_for_cross = bool(params.get("enforce_diff_pdb_for_cross", True))
hard_negative_topk = int(params.get("hard_negative_topk", 200))

parent_in_dir = params["masif_precomputation_dir"]
training_list = set(x.rstrip() for x in open(params["training_list"]).readlines())
testing_list = set(x.rstrip() for x in open(params["testing_list"]).readlines())
all_pairs = sorted(os.listdir(parent_in_dir))

print("Negative sampling config:")
print(f"  neg_ratio={neg_ratio}")
print(
    "  mix(cross,within,hard)=({:.3f}, {:.3f}, {:.3f})".format(
        mix[0], mix[1], mix[2]
    )
)
print(f"  enforce_diff_pdb_for_cross={enforce_diff_pdb_for_cross}")
print(f"  hard_negative_topk={hard_negative_topk}")

# Phase A: collect positives and per-pair negative candidates.
records = []
for count, ppi_pair_id in enumerate(all_pairs):
    print(f"{count} / {len(all_pairs)}")
    if ppi_pair_id not in training_list and ppi_pair_id not in testing_list:
        continue

    in_dir = os.path.join(parent_in_dir, ppi_pair_id) + "/"
    fields = ppi_pair_id.split("_")
    if len(fields) < 3:
        continue
    pdb_id = fields[0]
    split = "test"
    if ppi_pair_id in training_list:
        split = "train" if np.random.random() <= params["range_val_samples"] else "val"

    try:
        labels = np.load(in_dir + "p1_sc_labels.npy")
        labels = np.median(labels[0], axis=1)
    except Exception as e:
        print("Could not open {}p1_sc_labels.npy: {}".format(in_dir, e))
        continue

    ply_fn1 = masif_opts["ply_file_template"].format(fields[0], fields[1])
    ply_fn2 = masif_opts["ply_file_template"].format(fields[0], fields[2])

    pos_labels = np.where(
        (labels < params["max_sc_filt"]) & (labels > params["min_sc_filt"])
    )[0]
    k_accept = int(params["pos_surf_accept_probability"] * len(pos_labels))
    if k_accept < 1:
        continue
    chosen = np.arange(len(pos_labels))
    np.random.shuffle(chosen)
    chosen = pos_labels[chosen[:k_accept]]

    v1 = pymesh.load_mesh(ply_fn1).vertices[chosen]
    v2 = pymesh.load_mesh(ply_fn2).vertices

    kdt = cKDTree(v2)
    d, r = kdt.query(v1)
    contact_points = np.where(d < params["pos_interface_cutoff"])[0]
    if len(contact_points) == 0:
        continue
    k1 = chosen[contact_points]
    k2 = r[contact_points]
    assert len(k1) == len(k2)

    kdt = cKDTree(v1)
    dneg, _ = kdt.query(v2)
    k_neg2 = np.where(dneg > params["pos_interface_cutoff"])[0]
    if len(k_neg2) == 0:
        print(f"No non-interface negatives for {ppi_pair_id}, skipping.")
        continue

    p1_rho = np.load(in_dir + "p1_rho_wrt_center.npy")
    p1_theta = np.load(in_dir + "p1_theta_wrt_center.npy")
    p1_input = np.load(in_dir + "p1_input_feat.npy")
    p1_mask = np.load(in_dir + "p1_mask.npy")

    p2_rho = np.load(in_dir + "p2_rho_wrt_center.npy")
    p2_theta = np.load(in_dir + "p2_theta_wrt_center.npy")
    p2_input = np.load(in_dir + "p2_input_feat.npy")
    p2_mask = np.load(in_dir + "p2_mask.npy")

    record = {
        "ppi_pair_id": ppi_pair_id,
        "pdb_id": pdb_id,
        "split": split,
        "binder_rho": p1_rho[k1],
        "binder_theta": p1_theta[k1],
        "binder_input": p1_input[k1],
        "binder_mask": p1_mask[k1],
        "pos_rho": p2_rho[k2],
        "pos_theta": p2_theta[k2],
        "pos_input": p2_input[k2],
        "pos_mask": p2_mask[k2],
        "pos_names": [f"{ppi_pair_id}_p1_{ii}" for ii in k1],
        "within_rho": p2_rho[k_neg2],
        "within_theta": p2_theta[k_neg2],
        "within_input": p2_input[k_neg2],
        "within_mask": p2_mask[k_neg2],
        "within_names": [f"{ppi_pair_id}_p2_{ii}" for ii in k_neg2],
        "within_dneg": dneg[k_neg2],
    }
    records.append(record)

if len(records) == 0:
    raise RuntimeError("No valid records found to build cache.")

def build_pool(source_records, split_group):
    pool_rho = []
    pool_theta = []
    pool_input = []
    pool_mask = []
    pool_names = []
    pool_ppi = []
    pool_pdb = []
    for rec in source_records:
        if split_group == "trainval" and rec["split"] == "test":
            continue
        if split_group == "test" and rec["split"] != "test":
            continue
        n = len(rec["within_rho"])
        if n == 0:
            continue
        pool_rho.append(rec["within_rho"])
        pool_theta.append(rec["within_theta"])
        pool_input.append(rec["within_input"])
        pool_mask.append(rec["within_mask"])
        pool_names.extend(rec["within_names"])
        pool_ppi.extend([rec["ppi_pair_id"]] * n)
        pool_pdb.extend([rec["pdb_id"]] * n)
    return {
        "rho": concat_or_empty(pool_rho, 2),
        "theta": concat_or_empty(pool_theta, 2),
        "input": concat_or_empty(pool_input, 3),
        "mask": concat_or_empty(pool_mask, 2),
        "names": np.asarray(pool_names),
        "ppi": np.asarray(pool_ppi),
        "pdb": np.asarray(pool_pdb),
    }


pools = {
    "trainval": build_pool(records, "trainval"),
    "test": build_pool(records, "test"),
}

# Precompute cross-valid candidate indices per record.
for rec in records:
    pool_key = "test" if rec["split"] == "test" else "trainval"
    pool = pools[pool_key]
    cross_all = np.arange(len(pool["rho"]))
    if len(cross_all) == 0:
        rec["cross_valid_idx"] = np.asarray([], dtype=int)
    else:
        mask = pool["ppi"] != rec["ppi_pair_id"]
        if enforce_diff_pdb_for_cross:
            mask = mask & (pool["pdb"] != rec["pdb_id"])
        rec["cross_valid_idx"] = cross_all[mask]
    rec["cross_pool_key"] = pool_key
    rec["within_all_idx"] = np.arange(len(rec["within_rho"]))
    order = np.argsort(rec["within_dneg"])
    rec["hard_idx"] = order[: min(hard_negative_topk, len(order))]


# Phase B: mixed negative sampling and cache assembly.
binder_rho_wrt_center = []
binder_theta_wrt_center = []
binder_input_feat = []
binder_mask = []
pos_rho_wrt_center = []
pos_theta_wrt_center = []
pos_input_feat = []
pos_mask = []
neg_rho_wrt_center = []
neg_theta_wrt_center = []
neg_input_feat = []
neg_mask = []
pos_names = []
neg_names = []

pos_training_idx = []
pos_val_idx = []
pos_test_idx = []
neg_training_idx = []
neg_val_idx = []
neg_test_idx = []

mix_counts = {"cross": 0, "within": 0, "hard": 0, "fallback": 0}
pos_idx_count = 0
neg_idx_count = 0

for rec in records:
    n_pos = len(rec["pos_rho"])
    if n_pos == 0:
        continue

    # Positives (1 per binder-contact pair).
    binder_rho_wrt_center.append(rec["binder_rho"])
    binder_theta_wrt_center.append(rec["binder_theta"])
    binder_input_feat.append(rec["binder_input"])
    binder_mask.append(rec["binder_mask"])
    pos_rho_wrt_center.append(rec["pos_rho"])
    pos_theta_wrt_center.append(rec["pos_theta"])
    pos_input_feat.append(rec["pos_input"])
    pos_mask.append(rec["pos_mask"])
    pos_names.extend(rec["pos_names"])

    pos_indices = np.arange(pos_idx_count, pos_idx_count + n_pos)
    if rec["split"] == "train":
        pos_training_idx.extend(pos_indices.tolist())
    elif rec["split"] == "val":
        pos_val_idx.extend(pos_indices.tolist())
    else:
        pos_test_idx.extend(pos_indices.tolist())
    pos_idx_count += n_pos

    pool = pools[rec["cross_pool_key"]]
    cross_valid_idx = rec["cross_valid_idx"]
    cross_all_idx = np.arange(len(pool["rho"]))
    within_all_idx = rec["within_all_idx"]
    hard_idx = rec["hard_idx"]

    for _ in range(n_pos):
        bucket_counts = split_counts(neg_ratio, mix)
        sampled_neg_rows = []
        sampled_neg_names = []

        # Cross-complex bucket.
        n_cross = bucket_counts[0]
        cross_pick = sample_from_index_array(cross_valid_idx, n_cross)
        if len(cross_pick) > 0:
            sampled_neg_rows.append(
                (
                    pool["rho"][cross_pick],
                    pool["theta"][cross_pick],
                    pool["input"][cross_pick],
                    pool["mask"][cross_pick],
                )
            )
            sampled_neg_names.extend(pool["names"][cross_pick].tolist())
            mix_counts["cross"] += len(cross_pick)

        # Within-complex random bucket.
        n_within = bucket_counts[1]
        within_pick = sample_from_index_array(within_all_idx, n_within)
        if len(within_pick) > 0:
            sampled_neg_rows.append(
                (
                    rec["within_rho"][within_pick],
                    rec["within_theta"][within_pick],
                    rec["within_input"][within_pick],
                    rec["within_mask"][within_pick],
                )
            )
            sampled_neg_names.extend(np.asarray(rec["within_names"])[within_pick].tolist())
            mix_counts["within"] += len(within_pick)

        # Hard negatives (near-interface within-complex).
        n_hard = bucket_counts[2]
        hard_pick = sample_from_index_array(hard_idx, n_hard)
        if len(hard_pick) > 0:
            sampled_neg_rows.append(
                (
                    rec["within_rho"][hard_pick],
                    rec["within_theta"][hard_pick],
                    rec["within_input"][hard_pick],
                    rec["within_mask"][hard_pick],
                )
            )
            sampled_neg_names.extend(np.asarray(rec["within_names"])[hard_pick].tolist())
            mix_counts["hard"] += len(hard_pick)

        n_selected = len(sampled_neg_names)
        if n_selected < neg_ratio:
            n_missing = neg_ratio - n_selected
            fallback_source = "cross_valid"
            fallback_pick = sample_from_index_array(cross_valid_idx, n_missing)
            if len(fallback_pick) == 0:
                fallback_source = "within"
                fallback_pick = sample_from_index_array(within_all_idx, n_missing)
            if len(fallback_pick) == 0:
                fallback_source = "cross_all"
                fallback_pick = sample_from_index_array(cross_all_idx, n_missing)
            if len(fallback_pick) == 0:
                raise RuntimeError(
                    f"Could not sample fallback negatives for {rec['ppi_pair_id']}"
                )
            if fallback_source in ["cross_valid", "cross_all"]:
                sampled_neg_rows.append(
                    (
                        pool["rho"][fallback_pick],
                        pool["theta"][fallback_pick],
                        pool["input"][fallback_pick],
                        pool["mask"][fallback_pick],
                    )
                )
                sampled_neg_names.extend(pool["names"][fallback_pick].tolist())
            else:
                sampled_neg_rows.append(
                    (
                        rec["within_rho"][fallback_pick],
                        rec["within_theta"][fallback_pick],
                        rec["within_input"][fallback_pick],
                        rec["within_mask"][fallback_pick],
                    )
                )
                sampled_neg_names.extend(
                    np.asarray(rec["within_names"])[fallback_pick].tolist()
                )
            mix_counts["fallback"] += len(fallback_pick)

        # Materialize sampled negatives.
        for rho_chunk, theta_chunk, input_chunk, mask_chunk in sampled_neg_rows:
            neg_rho_wrt_center.append(rho_chunk)
            neg_theta_wrt_center.append(theta_chunk)
            neg_input_feat.append(input_chunk)
            neg_mask.append(mask_chunk)
        neg_names.extend(sampled_neg_names)

        n_new_neg = len(sampled_neg_names)
        neg_indices = np.arange(neg_idx_count, neg_idx_count + n_new_neg)
        if rec["split"] == "train":
            neg_training_idx.extend(neg_indices.tolist())
        elif rec["split"] == "val":
            neg_val_idx.extend(neg_indices.tolist())
        else:
            neg_test_idx.extend(neg_indices.tolist())
        neg_idx_count += n_new_neg

if not os.path.exists(params["cache_dir"]):
    os.makedirs(params["cache_dir"])

binder_rho_wrt_center = concat_or_empty(binder_rho_wrt_center, 2)
binder_theta_wrt_center = concat_or_empty(binder_theta_wrt_center, 2)
binder_input_feat = concat_or_empty(binder_input_feat, 3)
binder_mask = concat_or_empty(binder_mask, 2)
pos_rho_wrt_center = concat_or_empty(pos_rho_wrt_center, 2)
pos_theta_wrt_center = concat_or_empty(pos_theta_wrt_center, 2)
pos_input_feat = concat_or_empty(pos_input_feat, 3)
pos_mask = concat_or_empty(pos_mask, 2)
neg_rho_wrt_center = concat_or_empty(neg_rho_wrt_center, 2)
neg_theta_wrt_center = concat_or_empty(neg_theta_wrt_center, 2)
neg_input_feat = concat_or_empty(neg_input_feat, 3)
neg_mask = concat_or_empty(neg_mask, 2)

print(f"Positives before NaN filtering: {len(binder_input_feat)}")
print(f"Negatives before NaN filtering: {len(neg_input_feat)}")

pos_not_nan = ~np.isnan(binder_input_feat).any(axis=(1, 2))
pos_not_nan = pos_not_nan & (~np.isnan(pos_input_feat).any(axis=(1, 2)))
neg_not_nan = ~np.isnan(neg_input_feat).any(axis=(1, 2))

pos_names = [x for i, x in enumerate(pos_names) if pos_not_nan[i]]
neg_names = [x for i, x in enumerate(neg_names) if neg_not_nan[i]]
binder_rho_wrt_center = binder_rho_wrt_center[pos_not_nan, :]
binder_theta_wrt_center = binder_theta_wrt_center[pos_not_nan, :]
binder_input_feat = binder_input_feat[pos_not_nan, ...]
binder_mask = binder_mask[pos_not_nan, :]
pos_rho_wrt_center = pos_rho_wrt_center[pos_not_nan, :]
pos_theta_wrt_center = pos_theta_wrt_center[pos_not_nan, :]
pos_input_feat = pos_input_feat[pos_not_nan, ...]
pos_mask = pos_mask[pos_not_nan, :]
neg_rho_wrt_center = neg_rho_wrt_center[neg_not_nan, :]
neg_theta_wrt_center = neg_theta_wrt_center[neg_not_nan, :]
neg_input_feat = neg_input_feat[neg_not_nan, ...]
neg_mask = neg_mask[neg_not_nan, :]

pos_training_idx = remap_indices(pos_training_idx, pos_not_nan)
pos_val_idx = remap_indices(pos_val_idx, pos_not_nan)
pos_test_idx = remap_indices(pos_test_idx, pos_not_nan)
neg_training_idx = remap_indices(neg_training_idx, neg_not_nan)
neg_val_idx = remap_indices(neg_val_idx, neg_not_nan)
neg_test_idx = remap_indices(neg_test_idx, neg_not_nan)

assert len(set(pos_training_idx) & set(pos_val_idx)) == 0
assert len(set(pos_training_idx) & set(pos_test_idx)) == 0
assert len(set(pos_val_idx) & set(pos_test_idx)) == 0
assert len(set(neg_training_idx) & set(neg_val_idx)) == 0
assert len(set(neg_training_idx) & set(neg_test_idx)) == 0
assert len(set(neg_val_idx) & set(neg_test_idx)) == 0

if len(pos_training_idx) == 0 or len(neg_training_idx) == 0:
    raise RuntimeError("Insufficient train samples after filtering.")

print(f"Positives after NaN filtering: {len(binder_input_feat)}")
print(f"Negatives after NaN filtering: {len(neg_input_feat)}")
print(
    "Sampled negatives by bucket: cross={}, within={}, hard={}, fallback={}".format(
        mix_counts["cross"],
        mix_counts["within"],
        mix_counts["hard"],
        mix_counts["fallback"],
    )
)

np.save(params["cache_dir"] + "/pos_names.npy", pos_names)
np.save(params["cache_dir"] + "/neg_names.npy", neg_names)
np.save(params["cache_dir"] + "/binder_rho_wrt_center.npy", binder_rho_wrt_center)
np.save(params["cache_dir"] + "/binder_theta_wrt_center.npy", binder_theta_wrt_center)
np.save(params["cache_dir"] + "/binder_input_feat.npy", binder_input_feat)
np.save(params["cache_dir"] + "/binder_mask.npy", binder_mask)
np.save(params["cache_dir"] + "/pos_training_idx.npy", np.asarray(pos_training_idx, dtype=int))
np.save(params["cache_dir"] + "/pos_val_idx.npy", np.asarray(pos_val_idx, dtype=int))
np.save(params["cache_dir"] + "/pos_test_idx.npy", np.asarray(pos_test_idx, dtype=int))
np.save(params["cache_dir"] + "/pos_rho_wrt_center.npy", pos_rho_wrt_center)
np.save(params["cache_dir"] + "/pos_theta_wrt_center.npy", pos_theta_wrt_center)
np.save(params["cache_dir"] + "/pos_input_feat.npy", pos_input_feat)
np.save(params["cache_dir"] + "/pos_mask.npy", pos_mask)
np.save(params["cache_dir"] + "/neg_training_idx.npy", np.asarray(neg_training_idx, dtype=int))
np.save(params["cache_dir"] + "/neg_val_idx.npy", np.asarray(neg_val_idx, dtype=int))
np.save(params["cache_dir"] + "/neg_test_idx.npy", np.asarray(neg_test_idx, dtype=int))
np.save(params["cache_dir"] + "/neg_rho_wrt_center.npy", neg_rho_wrt_center)
np.save(params["cache_dir"] + "/neg_theta_wrt_center.npy", neg_theta_wrt_center)
np.save(params["cache_dir"] + "/neg_input_feat.npy", neg_input_feat)
np.save(params["cache_dir"] + "/neg_mask.npy", neg_mask)
