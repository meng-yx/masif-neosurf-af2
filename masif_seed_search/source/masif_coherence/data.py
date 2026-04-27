"""On-disk data loading for MaSIF-Coherence v1.

Consumes descriptor arrays and a mandatory prebuilt vertex cache:

    <descriptors_dir>/<pdbid>_<chains1>_<chains2>/
        p1_desc_straight.npy   # (V_A, 80)
        p1_desc_flipped.npy
        p2_desc_straight.npy   # (V_B, 80)
        p2_desc_flipped.npy

    <cache_dir>/manifest.json
    <cache_dir>/vertices/<one-file-per-pdb-chain>.npy

Positives are cognate pairs from the training/testing lists (one ppi_pair_id per
line). Negatives for a positive (p1_A, p2_A) are constructed by pairing p1_A
with p2_B sampled from a different ppi_pair_id.
"""

from __future__ import annotations

import json
import os
import random
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
from torch.utils.data import Dataset


# ---------------------------------------------------------------------------
# Single-protein loading.
# ---------------------------------------------------------------------------


@dataclass
class ProteinSurface:
    """Minimal v1 view of a protein surface: coordinates + MaSIF descriptors."""

    pdb_id: str
    chain: str
    side: str  # "p1" or "p2"
    desc_type: str  # "straight" or "flipped"
    X: torch.Tensor  # (V, 3) float32
    D: torch.Tensor  # (V, 80) float32


def _load_cache_entries(cache_dir: str) -> Dict[str, Dict[str, object]]:
    manifest_path = os.path.join(cache_dir, "manifest.json")
    if not os.path.exists(manifest_path):
        raise FileNotFoundError(
            f"Missing cache manifest at {manifest_path}. Run build_cache.py first."
        )
    with open(manifest_path, "r") as f:
        manifest = json.load(f)
    entries = manifest.get("entries")
    if not isinstance(entries, dict):
        raise RuntimeError(f"Malformed cache manifest at {manifest_path}.")
    return entries


def _load_cached_vertices(
    cache_dir: str, cache_entries: Dict[str, Dict[str, object]], chain_key: str
) -> np.ndarray:
    entry = cache_entries.get(chain_key)
    if entry is None:
        raise FileNotFoundError(
            f"No cache entry for chain {chain_key!r}. Re-run build_cache.py."
        )
    source_ply = str(entry["source_ply"])
    cached_rel = str(entry["cache_relpath"])
    cached_path = os.path.join(cache_dir, cached_rel)
    if not os.path.exists(cached_path):
        raise FileNotFoundError(
            f"Missing cached vertices file {cached_path}. Re-run build_cache.py."
        )
    if not os.path.exists(source_ply):
        raise FileNotFoundError(
            f"Source PLY does not exist anymore: {source_ply}. Re-run build_cache.py."
        )
    expected_mtime = float(entry["ply_mtime"])
    current_mtime = os.path.getmtime(source_ply)
    if abs(current_mtime - expected_mtime) > 1e-6:
        raise RuntimeError(
            f"Stale cache for {chain_key}. Source PLY changed: {source_ply}. "
            "Re-run build_cache.py."
        )
    verts = np.load(cached_path).astype(np.float32)
    if verts.ndim != 2 or verts.shape[1] != 3:
        raise RuntimeError(f"Unexpected cached vertex shape {verts.shape} at {cached_path}.")
    return verts


def load_protein(
    descriptors_dir: str,
    cache_dir: str,
    cache_entries: Dict[str, Dict[str, object]],
    ppi_pair_id: str,
    side: str,
    desc_type: str,
) -> ProteinSurface:
    """Load one side of a ppi_pair_id.

    Args:
        descriptors_dir: directory containing one sub-folder per ppi_pair_id
            with the `p{1,2}_desc_{straight,flipped}.npy` files.
        cache_dir: directory containing cached vertices and `manifest.json`.
        cache_entries: parsed map from `manifest.json["entries"]`.
        ppi_pair_id: e.g. `1A0G_A_B`.
        side: "p1" or "p2".
        desc_type: "straight" or "flipped".
    """
    if side not in ("p1", "p2"):
        raise ValueError(f"side must be p1 or p2, got {side!r}")
    if desc_type not in ("straight", "flipped"):
        raise ValueError(f"desc_type must be straight or flipped, got {desc_type!r}")

    fields = ppi_pair_id.split("_")
    if len(fields) < 3:
        raise ValueError(f"ppi_pair_id {ppi_pair_id!r} must be PDBID_CHAINS1_CHAINS2")
    pdb_id = fields[0]
    chain = fields[1] if side == "p1" else fields[2]

    pair_dir = os.path.join(descriptors_dir, ppi_pair_id)
    desc_path = os.path.join(pair_dir, f"{side}_desc_{desc_type}.npy")
    if not os.path.exists(desc_path):
        raise FileNotFoundError(desc_path)
    D = np.load(desc_path).astype(np.float32)
    if D.ndim != 2 or D.shape[1] != 80:
        raise RuntimeError(f"Unexpected descriptor shape {D.shape} at {desc_path}.")

    chain_key = f"{pdb_id}_{chain}"
    X = _load_cached_vertices(cache_dir, cache_entries, chain_key)

    # The .ply mesh vertex order is the canonical order MaSIF uses when producing
    # descriptors, so `X` and `D` should already be aligned row-for-row.
    if X.shape[0] != D.shape[0]:
        raise RuntimeError(
            f"Mesh has {X.shape[0]} vertices but descriptor has {D.shape[0]} rows"
            f" for {ppi_pair_id} side {side}."
        )

    return ProteinSurface(
        pdb_id=pdb_id,
        chain=chain,
        side=side,
        desc_type=desc_type,
        X=torch.from_numpy(X),
        D=torch.from_numpy(D),
    )


# ---------------------------------------------------------------------------
# Pair dataset.
# ---------------------------------------------------------------------------


def read_pair_list(path: str) -> List[str]:
    """Read a list of ppi_pair_ids, one per line, ignoring blanks and comments."""
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f.readlines()]
    return [ln for ln in lines if ln and not ln.startswith("#")]


class CoherenceDataset(Dataset):
    """Yields `(protein_A, protein_B, y)` tuples.

    Positives: `(p1, p2)` from the same `ppi_pair_id` (cognate interface).
    Negatives: `(p1 of pair_i, p2 of pair_j)` with `j != i`, drawn uniformly
    from the same list (per the v1 plan; no hard-negative mining in v1).

    Each positive is followed by `neg_per_pos` negatives in the dataset ordering,
    so with `neg_per_pos=1` and `len(pair_ids)=N` the dataset has `2*N` items.
    """

    def __init__(
        self,
        pair_ids: List[str],
        descriptors_dir: str,
        cache_dir: str,
        neg_per_pos: int = 1,
        seed: int = 0,
        skip_missing: bool = True,
    ):
        self.pair_ids = list(pair_ids)
        self.descriptors_dir = descriptors_dir
        self.cache_dir = cache_dir
        self.cache_entries = _load_cache_entries(cache_dir)
        self.neg_per_pos = int(neg_per_pos)
        self.skip_missing = skip_missing
        self._rng = random.Random(seed)
        self._surface_cache: Dict[Tuple[str, str, str], ProteinSurface] = {}

        # Pre-compute the flat index -> (pos_idx, is_positive, neg_source_idx) mapping.
        self._items: List[Tuple[int, int, Optional[int]]] = []
        for i in range(len(self.pair_ids)):
            self._items.append((i, 1, None))
            for _ in range(self.neg_per_pos):
                j = self._sample_neg_index(i)
                self._items.append((i, 0, j))

    def _sample_neg_index(self, exclude: int) -> int:
        if len(self.pair_ids) < 2:
            raise ValueError("Need at least 2 pair_ids to sample negatives.")
        while True:
            j = self._rng.randrange(len(self.pair_ids))
            if j != exclude:
                return j

    def __len__(self) -> int:
        return len(self._items)

    def _try_load(self, ppi_pair_id: str, side: str, desc_type: str):
        key = (ppi_pair_id, side, desc_type)
        if key in self._surface_cache:
            return self._surface_cache[key]
        try:
            protein = load_protein(
                self.descriptors_dir,
                self.cache_dir,
                self.cache_entries,
                ppi_pair_id,
                side,
                desc_type,
            )
            self._surface_cache[key] = protein
            return protein
        except (FileNotFoundError, RuntimeError) as exc:
            if self.skip_missing:
                return None
            raise exc

    def __getitem__(self, idx: int):
        pos_idx, label, neg_idx = self._items[idx]
        id_A = self.pair_ids[pos_idx]
        id_B = self.pair_ids[pos_idx if label == 1 else neg_idx]  # type: ignore[arg-type]
        # Convention from the v1 plan: A uses straight, B uses flipped.
        prot_A = self._try_load(id_A, "p1", "straight")
        prot_B = self._try_load(id_B, "p2", "flipped")
        if prot_A is None or prot_B is None:
            # Degrade gracefully: attempt the next index. This matches the MaSIF
            # training loops that `continue` past missing files.
            return self.__getitem__((idx + 1) % len(self))
        return prot_A, prot_B, torch.tensor(float(label), dtype=torch.float32)
