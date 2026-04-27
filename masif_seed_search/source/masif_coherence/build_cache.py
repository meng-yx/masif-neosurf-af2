"""Build mandatory vertex cache for MaSIF-Coherence training.

Usage:
    python build_cache.py <custom_params_module>
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib
import json
import os
import tempfile
import time
from typing import Dict, List, Set, Tuple

import numpy as np
from joblib import Parallel, delayed

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - fallback for minimal environments.
    tqdm = None

from config import get_default_params, merge_custom_params
from data import read_pair_list


def _parse_args():
    parser = argparse.ArgumentParser(description="Build cache_dir vertex cache for MaSIF-Coherence.")
    parser.add_argument(
        "custom_params_module",
        nargs="?",
        default="nn_models.v1.custom_params",
        help="Importable Python module containing custom_params dict.",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=os.cpu_count(),
        help="Number of parallel worker processes for PLY->NPY conversion.",
    )
    return parser.parse_args()


def _load_vertices_from_ply(ply_path: str) -> np.ndarray:
    import trimesh

    mesh = trimesh.load(ply_path, process=False)
    if not hasattr(mesh, "vertices"):
        raise RuntimeError(f"Loaded PLY at {ply_path} has no vertices attribute.")
    verts = np.asarray(mesh.vertices, dtype=np.float32)
    if verts.ndim != 2 or verts.shape[1] != 3:
        raise RuntimeError(f"Unexpected vertex array shape {verts.shape} at {ply_path}.")
    return verts


def _atomic_save_npy(path: str, arr: np.ndarray) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(prefix=".tmp_", suffix=".npy", dir=os.path.dirname(path))
    os.close(fd)
    try:
        np.save(tmp_path, arr)
        os.replace(tmp_path, path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def _parse_pair_id(ppi_pair_id: str) -> Tuple[str, str, str]:
    fields = ppi_pair_id.split("_")
    if len(fields) < 3:
        raise ValueError(f"Invalid ppi_pair_id {ppi_pair_id!r}, expected PDBID_CHAIN1_CHAIN2")
    return fields[0], fields[1], fields[2]


def _cache_filename(chain_key: str, ply_path: str) -> str:
    digest = hashlib.sha1(os.path.abspath(ply_path).encode("utf-8")).hexdigest()[:12]
    return f"{chain_key}__{digest}.npy"


def _resolve_cache_dir(params: Dict[str, object]) -> str:
    if "model_dir" not in params:
        raise KeyError("Missing required parameter 'model_dir'.")
    return os.path.join(str(params["model_dir"]), "cache")


def _process_chain(ply_chain_dir: str, vertices_dir: str, pdb_id: str, chain: str):
    chain_key = f"{pdb_id}_{chain}"
    ply_path = os.path.join(ply_chain_dir, f"{chain_key}.ply")
    if not os.path.exists(ply_path):
        return {
            "ok": False,
            "chain_key": chain_key,
            "source_ply": os.path.abspath(ply_path),
            "reason": "missing_ply",
            "error": f"Missing source PLY file: {ply_path}",
        }

    try:
        verts = _load_vertices_from_ply(ply_path)
        cache_name = _cache_filename(chain_key, ply_path)
        cache_path = os.path.join(vertices_dir, cache_name)
        _atomic_save_npy(cache_path, verts)

        entry = {
            "source_ply": os.path.abspath(ply_path),
            "cache_relpath": os.path.join("vertices", cache_name),
            "ply_mtime": os.path.getmtime(ply_path),
            "shape": list(verts.shape),
            "dtype": str(verts.dtype),
        }
        return {"ok": True, "chain_key": chain_key, "entry": entry}
    except Exception as exc:
        return {
            "ok": False,
            "chain_key": chain_key,
            "source_ply": os.path.abspath(ply_path),
            "reason": "cache_build_error",
            "error": repr(exc),
        }


@contextlib.contextmanager
def _tqdm_joblib(total: int):
    from joblib import parallel

    if tqdm is None:
        yield
        return

    pbar = tqdm(total=total, desc="Building vertex cache", unit="chain")
    old_callback = parallel.BatchCompletionCallBack

    class TqdmBatchCompletionCallback(old_callback):
        def __call__(self, *args, **kwargs):
            pbar.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield
    finally:
        parallel.BatchCompletionCallBack = old_callback
        pbar.close()


def main() -> None:
    args = _parse_args()
    params = get_default_params()
    custom_mod = importlib.import_module(args.custom_params_module, package=None)
    params = merge_custom_params(params, custom_mod.custom_params)

    for req_key in ("model_dir", "ply_chain_dir", "training_list", "testing_list"):
        if req_key not in params:
            raise KeyError(f"Missing required parameter {req_key!r}.")

    cache_dir = _resolve_cache_dir(params)
    vertices_dir = os.path.join(cache_dir, "vertices")
    os.makedirs(vertices_dir, exist_ok=True)

    pair_ids: List[str] = []
    pair_ids.extend(read_pair_list(params["training_list"]))
    pair_ids.extend(read_pair_list(params["testing_list"]))
    pair_ids = sorted(set(pair_ids))

    required_chains: Set[Tuple[str, str]] = set()
    for pair_id in pair_ids:
        pdb_id, chain_1, chain_2 = _parse_pair_id(pair_id)
        required_chains.add((pdb_id, chain_1))
        required_chains.add((pdb_id, chain_2))

    sorted_chains = sorted(required_chains)
    num_workers = max(1, int(args.num_workers))
    print(f"Building cache for {len(sorted_chains)} chains with {num_workers} workers")
    with _tqdm_joblib(total=len(sorted_chains)):
        records = Parallel(n_jobs=num_workers, prefer="processes")(
            delayed(_process_chain)(params["ply_chain_dir"], vertices_dir, pdb_id, chain)
            for pdb_id, chain in sorted_chains
        )
    entries: Dict[str, Dict[str, object]] = {}
    missing_chains: List[Dict[str, object]] = []
    for record in records:
        if record["ok"]:
            entries[record["chain_key"]] = record["entry"]  # type: ignore[index]
        else:
            missing_chains.append(
                {
                    "chain_key": record["chain_key"],
                    "source_ply": record["source_ply"],
                    "reason": record["reason"],
                    "error": record["error"],
                }
            )

    if missing_chains:
        print(f"WARNING: encountered {len(missing_chains)} missing/unusable chains.")
        for row in missing_chains[:10]:
            print(
                "  missing chain {chain_key} reason={reason} source={source_ply}".format(
                    **row
                )
            )
        if len(missing_chains) > 10:
            print(f"  ... and {len(missing_chains) - 10} more.")

    manifest = {
        "version": 1,
        "created_at_epoch_s": time.time(),
        "custom_params_module": args.custom_params_module,
        "pair_count": len(pair_ids),
        "chain_count": len(required_chains),
        "requested_chain_count": len(sorted_chains),
        "cached_chain_count": len(entries),
        "missing_chain_count": len(missing_chains),
        "missing_chains": missing_chains,
        "entries": entries,
    }
    manifest_path = os.path.join(cache_dir, "manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    print(f"Wrote vertex cache manifest to {manifest_path}")


if __name__ == "__main__":
    main()
