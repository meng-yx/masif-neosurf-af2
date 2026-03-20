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
        "--catalog-dir",
        default=None,
        help="Catalog output directory. Defaults to <cache_dir>/catalog.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional numpy random seed for deterministic split assignment.",
    )
    return parser.parse_args()


def build_record(params, ppi_pair_id):
    in_dir = os.path.join(params["masif_precomputation_dir"], ppi_pair_id)
    fields = ppi_pair_id.split("_")
    if len(fields) < 3:
        return None
    pdb_id = fields[0]

    if ppi_pair_id in params["_training_list"]:
        split = "train" if np.random.random() <= params["range_val_samples"] else "val"
    elif ppi_pair_id in params["_testing_list"]:
        split = "test"
    else:
        return None

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
    catalog_dir = args.catalog_dir or os.path.join(params["cache_dir"], "catalog")
    ensure_dir(catalog_dir)

    print("Building cache catalog in {}".format(catalog_dir))
    print("Scanning {} candidate directories".format(len(all_pairs)))
    records = []
    for count, ppi_pair_id in enumerate(all_pairs):
        if count % 100 == 0:
            print("{}/{}".format(count, len(all_pairs)))
        rec = build_record(params, ppi_pair_id)
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
