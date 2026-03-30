import argparse
import os
import time

import numpy as np
import pymesh
from scipy.spatial import cKDTree

from cache_shard_common import (
    ensure_dir,
    load_params,
    params_fingerprint,
    save_json,
    write_success_marker,
)
from default_config.masif_opts import masif_opts


def parse_args():
    parser = argparse.ArgumentParser(description="Build metadata catalog for sharded cache.")
    parser.add_argument(
        "custom_params_module",
        nargs="?",
        default=None,
        help="Python module path containing custom_params dictionary",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional numpy random seed for deterministic split assignment.",
    )
    return parser.parse_args()


def resolve_split(params, ppi_pair_id):
    if ppi_pair_id in params["_training_list"]:
        return "train" if np.random.random() <= params["range_val_samples"] else "val"
    elif ppi_pair_id in params["_testing_list"]:
        return "test"
    return None


def get_model_type(ppi_pair_id):
    # Example: 1ATP-PP_L_R -> PP
    try:
        return ppi_pair_id.split("-")[1].split("_")[0]
    except Exception:
        return None


def to_pp_complex_id(ppi_pair_id):
    model_type = get_model_type(ppi_pair_id)
    if model_type is None:
        return None
    return ppi_pair_id.replace("-{}_".format(model_type), "-PP_", 1)


def build_pp_record(params, ppi_pair_id, split):
    in_dir = os.path.join(params["masif_precomputation_dir"], ppi_pair_id)
    fields = ppi_pair_id.split("_")
    if len(fields) < 3:
        return None
    pdb_id = fields[0]

    try:
        labels = np.load(os.path.join(in_dir, "p1_sc_labels.npy"))
        labels = np.median(labels[0], axis=1)
    except Exception as exc:
        print("Could not open {}: {}".format(os.path.join(in_dir, "p1_sc_labels.npy"), exc))
        return None

    ply_fn1 = masif_opts["ply_file_template"].format(fields[0], fields[1])
    ply_fn2 = masif_opts["ply_file_template"].format(fields[0], fields[2])

    pos_labels = np.where(
        (labels < params["max_sc_filt"]) & (labels > params["min_sc_filt"])
    )[0]
    k_accept = int(params["pos_surf_accept_probability"] * len(pos_labels))
    if k_accept < 1:
        return None
    chosen = np.arange(len(pos_labels))
    np.random.shuffle(chosen)
    chosen = pos_labels[chosen[:k_accept]]

    v1 = pymesh.load_mesh(ply_fn1).vertices[chosen]
    v2 = pymesh.load_mesh(ply_fn2).vertices

    kdt = cKDTree(v2)
    d, r = kdt.query(v1)
    contact_points = np.where(d < params["pos_interface_cutoff"])[0]
    if len(contact_points) == 0:
        return None

    k1 = chosen[contact_points]
    k2 = r[contact_points]

    kdt = cKDTree(v1)
    dneg, _ = kdt.query(v2)
    k_neg2 = np.where(dneg > params["pos_interface_cutoff"])[0]
    if len(k_neg2) == 0:
        print("No non-interface negatives for {}, skipping.".format(ppi_pair_id))
        return None

    return {
        "ppi_pair_id": ppi_pair_id,
        "pdb_id": pdb_id,
        "split": split,
        "in_dir": in_dir,
        "k1": np.asarray(k1, dtype=int),
        "k2": np.asarray(k2, dtype=int),
        "k_neg2": np.asarray(k_neg2, dtype=int),
        "within_dneg": np.asarray(dneg[k_neg2], dtype=float),
    }


def build_mapped_record_from_pp(params, target_ppi_pair_id, split, pp_record):
    target_in_dir = os.path.join(params["masif_precomputation_dir"], target_ppi_pair_id)
    target_fields = target_ppi_pair_id.split("_")
    if len(target_fields) < 3:
        return None

    pp_fields = pp_record["ppi_pair_id"].split("_")
    if len(pp_fields) < 3:
        return None

    try:
        pp_ply1 = masif_opts["ply_file_template"].format(pp_fields[0], pp_fields[1])
        pp_ply2 = masif_opts["ply_file_template"].format(pp_fields[0], pp_fields[2])
        tgt_ply1 = masif_opts["ply_file_template"].format(target_fields[0], target_fields[1])
        tgt_ply2 = masif_opts["ply_file_template"].format(target_fields[0], target_fields[2])

        pp_v1 = pymesh.load_mesh(pp_ply1).vertices
        pp_v2 = pymesh.load_mesh(pp_ply2).vertices
        target_v1 = pymesh.load_mesh(tgt_ply1).vertices
        target_v2 = pymesh.load_mesh(tgt_ply2).vertices
    except Exception as exc:
        print("Could not load mapping meshes for {}: {}".format(target_ppi_pair_id, exc))
        return None

    if len(pp_record["k1"]) == 0 or len(pp_record["k2"]) == 0:
        return None

    pp_coords_p1 = pp_v1[np.asarray(pp_record["k1"], dtype=int)]
    pp_coords_p2 = pp_v2[np.asarray(pp_record["k2"], dtype=int)]

    kdt_tgt_v1 = cKDTree(target_v1)
    _, k1_map = kdt_tgt_v1.query(pp_coords_p1)
    k1_map = np.asarray(k1_map, dtype=int)

    kdt_tgt_v2 = cKDTree(target_v2)
    _, k2_map = kdt_tgt_v2.query(pp_coords_p2)
    k2_map = np.asarray(k2_map, dtype=int)

    if len(k1_map) == 0 or len(k2_map) == 0:
        return None

    v1_candidates = target_v1[k1_map]
    kdt_v1 = cKDTree(v1_candidates)
    dneg, _ = kdt_v1.query(target_v2)
    k_neg2 = np.where(dneg > params["pos_interface_cutoff"])[0]
    if len(k_neg2) == 0:
        print("No non-interface negatives for {}, skipping.".format(target_ppi_pair_id))
        return None

    return {
        "ppi_pair_id": target_ppi_pair_id,
        "pdb_id": target_fields[0],
        "split": split,
        "in_dir": target_in_dir,
        "k1": k1_map,
        "k2": k2_map,
        "k_neg2": np.asarray(k_neg2, dtype=int),
        "within_dneg": np.asarray(dneg[k_neg2], dtype=float),
    }


def build_record(params, ppi_pair_id, pp_record_cache):
    split = resolve_split(params, ppi_pair_id)
    if split is None:
        return None

    model_type = get_model_type(ppi_pair_id)
    if model_type == "PP":
        print("Building PP-native positives for {}".format(ppi_pair_id))
        return build_pp_record(params, ppi_pair_id, split)

    if model_type in ["PA", "AP"]:
        pp_id = to_pp_complex_id(ppi_pair_id)
        if pp_id not in pp_record_cache:
            pp_split = resolve_split(params, pp_id)
            if pp_split is None:
                print("Could not resolve split for PP source {}, skipping {}.".format(pp_id, ppi_pair_id))
                return None
            pp_record_cache[pp_id] = build_pp_record(params, pp_id, pp_split)
        pp_record = pp_record_cache.get(pp_id)
        if pp_record is None:
            print("Could not build PP source {}, skipping {}.".format(pp_id, ppi_pair_id))
            return None
        print("Building mapped-from-PP positives for {} via {}".format(ppi_pair_id, pp_id))
        return build_mapped_record_from_pp(params, ppi_pair_id, split, pp_record)

    # Leave AA (and any unknown model type) unchanged for now.
    return None


def build_pools(records):
    pools = {
        "trainval": {
            "source_record_idx": [],
            "local_within_idx": [],
        },
        "test": {
            "source_record_idx": [],
            "local_within_idx": [],
        },
    }
    for rec_idx, rec in enumerate(records):
        pool_key = "test" if rec["split"] == "test" else "trainval"
        for local_idx, _ in enumerate(rec["k_neg2"]):
            pools[pool_key]["source_record_idx"].append(rec_idx)
            pools[pool_key]["local_within_idx"].append(local_idx)
    for pool_key in pools:
        for field in ["source_record_idx", "local_within_idx"]:
            pools[pool_key][field] = np.asarray(pools[pool_key][field], dtype=np.int32)
    return pools


def collect_allowed_pairs(params):
    """
    Restrict catalog candidates to:
      1) PP complexes explicitly listed in training/testing lists.
      2) Existing non-PP variants (PA/AP/AA/...) whose PP source is listed.
    """
    listed_ids = params["_training_list"] | params["_testing_list"]
    pp_ids = sorted(
        ppi_pair_id for ppi_pair_id in listed_ids if get_model_type(ppi_pair_id) == "PP"
    )
    pp_id_set = set(pp_ids)

    precomp_dir = params["masif_precomputation_dir"]
    existing_dirs = sorted(os.listdir(precomp_dir))
    existing_set = set(existing_dirs)

    missing_pp = sorted(pp_id for pp_id in pp_ids if pp_id not in existing_set)
    if missing_pp:
        print(
            "Warning: {} PP IDs listed but missing in precomputation dir.".format(
                len(missing_pp)
            )
        )
        print("First 10 missing PP IDs: {}".format(missing_pp[:10]))

    available_pp = [pp_id for pp_id in pp_ids if pp_id in existing_set]
    derived_ids = [
        ppi_pair_id
        for ppi_pair_id in existing_dirs
        if get_model_type(ppi_pair_id) != "PP"
        and to_pp_complex_id(ppi_pair_id) in pp_id_set
    ]

    all_pairs = sorted(set(available_pp) | set(derived_ids))
    return all_pairs, available_pp, derived_ids


def main():
    args = parse_args()
    params = load_params(args.custom_params_module)
    if args.seed is not None:
        np.random.seed(args.seed)

    params["_training_list"] = set(
        x.rstrip() for x in open(params["training_list"]).readlines()
    )
    params["_testing_list"] = set(
        x.rstrip() for x in open(params["testing_list"]).readlines()
    )
    all_pairs, available_pp, derived_ids = collect_allowed_pairs(params)
    # Keep all cache artifacts under params["cache_dir"] (custom_params.py).
    catalog_dir = os.path.join(params["cache_dir"], "catalog")
    ensure_dir(catalog_dir)

    print("Building cache catalog in {}".format(catalog_dir))
    print(
        "Restricted candidate set from train/test lists: "
        "{} PP + {} mapped variants = {} total".format(
            len(available_pp), len(derived_ids), len(all_pairs)
        )
    )
    if len(all_pairs) == 0:
        raise RuntimeError(
            "No candidate pairs found after list-based filtering. "
            "Check training/testing lists and masif_precomputation_dir."
        )
    records = []
    pp_record_cache = {}
    for count, ppi_pair_id in enumerate(all_pairs):
        if count % 100 == 0:
            print("{}/{}".format(count, len(all_pairs)))
        rec = build_record(params, ppi_pair_id, pp_record_cache)
        if rec is not None:
            records.append(rec)
    if len(records) == 0:
        raise RuntimeError("No valid records found to build catalog.")

    pools = build_pools(records)
    records_path = os.path.join(catalog_dir, "records.npy")
    np.save(records_path, np.asarray(records, dtype=object), allow_pickle=True)
    for pool_key, pool in pools.items():
        np.save(
            os.path.join(catalog_dir, "{}_source_record_idx.npy".format(pool_key)),
            pool["source_record_idx"],
        )
        np.save(
            os.path.join(catalog_dir, "{}_local_within_idx.npy".format(pool_key)),
            pool["local_within_idx"],
        )

    fingerprint_payload = {
        "training_list_path": params["training_list"],
        "testing_list_path": params["testing_list"],
        "neg_ratio": int(params.get("neg_ratio", 1)),
        "neg_mix_cross_complex": float(params.get("neg_mix_cross_complex", 0.0)),
        "neg_mix_within_complex": float(params.get("neg_mix_within_complex", 1.0)),
        "neg_mix_hard": float(params.get("neg_mix_hard", 0.0)),
        "enforce_diff_pdb_for_cross": bool(
            params.get("enforce_diff_pdb_for_cross", True)
        ),
        "hard_negative_topk": int(params.get("hard_negative_topk", 200)),
        "range_val_samples": float(params["range_val_samples"]),
        "max_sc_filt": float(params["max_sc_filt"]),
        "min_sc_filt": float(params["min_sc_filt"]),
        "pos_surf_accept_probability": float(params["pos_surf_accept_probability"]),
        "pos_interface_cutoff": float(params["pos_interface_cutoff"]),
        "masif_precomputation_dir": params["masif_precomputation_dir"],
        "cache_dir": params["cache_dir"],
        "seed": args.seed,
    }
    fingerprint = params_fingerprint(fingerprint_payload)
    manifest = {
        "version": 1,
        "created_at_epoch_s": time.time(),
        "fingerprint": fingerprint,
        "records_count": len(records),
        "trainval_pool_count": int(len(pools["trainval"]["source_record_idx"])),
        "test_pool_count": int(len(pools["test"]["source_record_idx"])),
        "fingerprint_payload": fingerprint_payload,
    }
    save_json(os.path.join(catalog_dir, "manifest.json"), manifest)
    write_success_marker(catalog_dir)
    print("Catalog complete with {} records.".format(len(records)))


if __name__ == "__main__":
    main()
