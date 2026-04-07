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


def _bump(stats, key):
    if stats is not None:
        stats[key] = stats.get(key, 0) + 1


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
    # Expected ID format: <pdb>-<chainL>-<chainR>-<TYPE>_<chain1>_<chain2>
    # Example: 1ATP-LA-LB-PP_L_R -> PP
    model_token = ppi_pair_id.split("_")[0]
    model_fields = model_token.split("-")
    if len(model_fields) >= 4:
        return model_fields[3]
    return None


def to_pp_complex_id(ppi_pair_id):
    model_type = get_model_type(ppi_pair_id)
    if model_type is None:
        return ppi_pair_id
    return ppi_pair_id.replace("-{}_".format(model_type), "-PP_", 1)


def build_pp_record(params, ppi_pair_id, split, stats=None):
    in_dir = os.path.join(params["masif_precomputation_dir"], ppi_pair_id)
    fields = ppi_pair_id.split("_")
    if len(fields) < 3:
        _bump(stats, "pp_bad_id_format")
        return None
    pdb_id = fields[0]

    try:
        labels = np.load(os.path.join(in_dir, "p1_sc_labels.npy"))
        labels = np.median(labels[0], axis=1)
    except Exception as exc:
        print("Could not open {}: {}".format(os.path.join(in_dir, "p1_sc_labels.npy"), exc))
        _bump(stats, "pp_missing_or_bad_labels")
        return None

    ply_fn1 = masif_opts["ply_file_template"].format(fields[0], fields[1])
    ply_fn2 = masif_opts["ply_file_template"].format(fields[0], fields[2])

    pos_labels = np.where(
        (labels < params["max_sc_filt"]) & (labels > params["min_sc_filt"])
    )[0]
    k_accept = int(params["pos_surf_accept_probability"] * len(pos_labels))
    if k_accept < 1:
        _bump(stats, "pp_no_positive_labels_after_filter")
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
        _bump(stats, "pp_no_interface_contacts")
        return None

    k1 = chosen[contact_points]
    k2 = r[contact_points]

    kdt = cKDTree(v1)
    dneg, _ = kdt.query(v2)
    k_neg2 = np.where(dneg > params["pos_interface_cutoff"])[0]
    if len(k_neg2) == 0:
        print("No non-interface negatives for {}, skipping.".format(ppi_pair_id))
        _bump(stats, "pp_no_noninterface_negatives")
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


def build_mapped_record_from_pp(params, target_ppi_pair_id, split, pp_record, stats=None):
    target_in_dir = os.path.join(params["masif_precomputation_dir"], target_ppi_pair_id)
    target_fields = target_ppi_pair_id.split("_")
    if len(target_fields) < 3:
        _bump(stats, "map_bad_target_id_format")
        return None

    pp_fields = pp_record["ppi_pair_id"].split("_")
    if len(pp_fields) < 3:
        _bump(stats, "map_bad_source_pp_id_format")
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
        _bump(stats, "map_mesh_load_error")
        return None

    if len(pp_record["k1"]) == 0 or len(pp_record["k2"]) == 0:
        _bump(stats, "map_empty_pp_indices")
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
        _bump(stats, "map_empty_mapped_indices")
        return None

    v1_candidates = target_v1[k1_map]
    kdt_v1 = cKDTree(v1_candidates)
    dneg, _ = kdt_v1.query(target_v2)
    k_neg2 = np.where(dneg > params["pos_interface_cutoff"])[0]
    if len(k_neg2) == 0:
        print("No non-interface negatives for {}, skipping.".format(target_ppi_pair_id))
        _bump(stats, "map_no_noninterface_negatives")
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


def build_record(params, ppi_pair_id, pp_record_cache, stats=None):
    split = resolve_split(params, ppi_pair_id)
    if split is None:
        _bump(stats, "split_not_listed")
        return None

    model_type = get_model_type(ppi_pair_id)
    if model_type == "PP":
        print("Building PP-native positives for {}".format(ppi_pair_id))
        rec = build_pp_record(params, ppi_pair_id, split, stats=stats)
        if rec is not None:
            _bump(stats, "built_pp")
        return rec

    if model_type in ["PA", "AP"]:
        pp_id = to_pp_complex_id(ppi_pair_id)
        if pp_id not in pp_record_cache:
            pp_split = resolve_split(params, pp_id)
            if pp_split is None:
                # PP source may be intentionally absent from train/test lists.
                # For mapped AP/PA entries listed in train/test, inherit the target split.
                pp_split = split
                print(
                    "PP source {} not listed in split files; using target split {} for {}.".format(
                        pp_id, split, ppi_pair_id
                    )
                )
            pp_record_cache[pp_id] = build_pp_record(params, pp_id, pp_split, stats=stats)
        pp_record = pp_record_cache.get(pp_id)
        if pp_record is None:
            print("Could not build PP source {}, skipping {}.".format(pp_id, ppi_pair_id))
            _bump(stats, "map_missing_pp_source_record")
            return None
        print("Building mapped-from-PP positives for {} via {}".format(ppi_pair_id, pp_id))
        rec = build_mapped_record_from_pp(
            params, ppi_pair_id, split, pp_record, stats=stats
        )
        if rec is not None:
            _bump(stats, "built_mapped")
        return rec

    # Leave AA (and any unknown model type) unchanged for now.
    _bump(stats, "unsupported_model_type")
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
    all_pairs = sorted(os.listdir(params["masif_precomputation_dir"]))
    # Keep all cache artifacts under params["cache_dir"] (custom_params.py).
    catalog_dir = os.path.join(params["cache_dir"], "catalog")
    ensure_dir(catalog_dir)

    print("Building cache catalog in {}".format(catalog_dir))
    print("Scanning {} candidate directories".format(len(all_pairs)))
    records = []
    stats = {}
    pp_record_cache = {}
    for count, ppi_pair_id in enumerate(all_pairs):
        if count % 100 == 0:
            print("{}/{}".format(count, len(all_pairs)))
        rec = build_record(params, ppi_pair_id, pp_record_cache, stats=stats)
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
    total_candidates = len(all_pairs)
    failed_total = total_candidates - len(records)
    print("Catalog summary: {} scanned, {} built, {} failed.".format(total_candidates, len(records), failed_total))
    if failed_total > 0:
        print("Failure reason counts:")
        ignore_keys = {"built_pp", "built_mapped"}
        reason_counts = [
            (k, int(v))
            for k, v in stats.items()
            if k not in ignore_keys and int(v) > 0
        ]
        for key, value in sorted(reason_counts, key=lambda x: (-x[1], x[0])):
            print("  - {}: {}".format(key, value))


if __name__ == "__main__":
    main()
