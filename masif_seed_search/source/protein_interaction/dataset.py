#!/usr/bin/env python3
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np

from protein_interaction.utils import resolve_path


@dataclass
class PairRecord:
    query_id: str
    matched_id: str
    label: int
    split: str


def parse_pair_id_chain(ppi_id: str, chain_index: int) -> str:
    fields = ppi_id.split("_")
    if len(fields) < chain_index + 1:
        raise ValueError(f"Cannot parse chain index {chain_index} from {ppi_id}")
    pdb_id = fields[0]
    chain_id = fields[chain_index]
    return f"{pdb_id}_{chain_id}"


def read_pairs_csv(path_like: str) -> List[PairRecord]:
    path = resolve_path(path_like)
    rows: List[PairRecord] = []
    with open(path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                PairRecord(
                    query_id=row["query_id"],
                    matched_id=row["matched_id"],
                    label=int(row["label"]),
                    split=row.get("split", ""),
                )
            )
    return rows


def load_vertices_from_ply(surfaces_dir: Path, pdb_chain: str) -> np.ndarray:
    ply_path = surfaces_dir / f"{pdb_chain}.ply"
    if not ply_path.exists():
        raise FileNotFoundError(str(ply_path))
    with open(ply_path, "r") as f:
        lines = f.readlines()

    n_verts = 0
    header_end = 0
    for i, line in enumerate(lines):
        if line.startswith("element vertex"):
            n_verts = int(line.split()[2])
        if line.strip() == "end_header":
            header_end = i + 1
            break
    if n_verts <= 0:
        raise ValueError(f"Invalid ply vertices in {ply_path}")

    verts = []
    for i in range(header_end, header_end + n_verts):
        fields = lines[i].strip().split()
        verts.append([float(fields[0]), float(fields[1]), float(fields[2])])
    return np.asarray(verts, dtype=np.float32)


def load_iface_scores(iface_dir: Path, pdb_chain: str) -> np.ndarray:
    pred_path = iface_dir / f"pred_{pdb_chain}.npy"
    if not pred_path.exists():
        return np.array([], dtype=np.float32)
    iface = np.load(pred_path)
    iface = np.squeeze(iface)
    return iface.astype(np.float32)


class ProteinPairDataset:
    def __init__(self, cfg: Dict, pair_records: Iterable[PairRecord], seed: int = 42):
        self.cfg = cfg
        self.records = list(pair_records)
        self.rng = np.random.RandomState(seed)

        data_cfg = cfg["data"]
        self.descriptors_dir = resolve_path(data_cfg["descriptors_dir"])
        self.surfaces_dir = resolve_path(data_cfg["surfaces_dir"])
        self.iface_dir = resolve_path(data_cfg["iface_dir"])
        self.query_desc_file = data_cfg["query_desc_file"]
        self.matched_desc_file = data_cfg["matched_desc_file"]
        self.query_chain_index = int(data_cfg.get("query_chain_index", 1))
        self.matched_chain_index = int(data_cfg.get("matched_chain_index", 2))
        self.use_query_interface_subset = bool(data_cfg.get("use_query_interface_subset", True))
        self.iface_cutoff = float(data_cfg.get("iface_cutoff", 0.7))
        self.query_iface_fallback_topk = int(data_cfg.get("query_iface_fallback_topk", 256))
        self.max_query_points = int(data_cfg["max_query_points"])
        self.max_matched_points = int(data_cfg["max_matched_points"])

    def __len__(self):
        return len(self.records)

    def _subsample_idx(self, n: int, max_points: int) -> np.ndarray:
        if n <= max_points:
            return np.arange(n, dtype=int)
        return self.rng.choice(n, size=max_points, replace=False)

    def _query_indices(self, query_id: str, n_vertices: int) -> np.ndarray:
        if not self.use_query_interface_subset:
            return np.arange(n_vertices, dtype=int)

        pdb_chain = parse_pair_id_chain(query_id, self.query_chain_index)
        iface = load_iface_scores(self.iface_dir, pdb_chain)
        if iface.size == 0 or len(iface) != n_vertices:
            return np.arange(n_vertices, dtype=int)

        idx = np.where(iface > self.iface_cutoff)[0]
        if len(idx) > 0:
            return idx

        k = min(self.query_iface_fallback_topk, n_vertices)
        if k == n_vertices:
            return np.arange(n_vertices, dtype=int)
        top_idx = np.argsort(iface)[-k:]
        return np.sort(top_idx)

    def _load_query(self, query_id: str) -> Tuple[np.ndarray, np.ndarray]:
        desc_path = self.descriptors_dir / query_id / self.query_desc_file
        query_desc = np.load(desc_path).astype(np.float32)
        query_chain = parse_pair_id_chain(query_id, self.query_chain_index)
        try:
            query_xyz = load_vertices_from_ply(self.surfaces_dir, query_chain)
        except FileNotFoundError:
            # Fallback for environments where benchmark surfaces are not present.
            query_xyz = np.zeros((len(query_desc), 3), dtype=np.float32)

        n = min(len(query_desc), len(query_xyz))
        query_desc = query_desc[:n]
        query_xyz = query_xyz[:n]

        idx = self._query_indices(query_id, n)
        query_desc = query_desc[idx]
        query_xyz = query_xyz[idx]

        keep_idx = self._subsample_idx(len(query_desc), self.max_query_points)
        query_desc = query_desc[keep_idx]
        query_xyz = query_xyz[keep_idx]
        return query_desc, query_xyz

    def _load_matched(self, matched_id: str) -> Tuple[np.ndarray, np.ndarray]:
        desc_path = self.descriptors_dir / matched_id / self.matched_desc_file
        matched_desc = np.load(desc_path).astype(np.float32)
        matched_chain = parse_pair_id_chain(matched_id, self.matched_chain_index)
        try:
            matched_xyz = load_vertices_from_ply(self.surfaces_dir, matched_chain)
        except FileNotFoundError:
            matched_xyz = np.zeros((len(matched_desc), 3), dtype=np.float32)

        n = min(len(matched_desc), len(matched_xyz))
        matched_desc = matched_desc[:n]
        matched_xyz = matched_xyz[:n]

        if len(matched_desc) > self.max_matched_points:
            idx = self._subsample_idx(len(matched_desc), self.max_matched_points)
            matched_desc = matched_desc[idx]
            matched_xyz = matched_xyz[idx]

        return matched_desc, matched_xyz

    def __getitem__(self, idx: int):
        row = self.records[idx]
        q_desc, q_xyz = self._load_query(row.query_id)
        m_desc, m_xyz = self._load_matched(row.matched_id)
        return {
            "query_desc": q_desc,
            "query_xyz": q_xyz,
            "matched_desc": m_desc,
            "matched_xyz": m_xyz,
            "label": np.int64(row.label),
            "query_id": row.query_id,
            "matched_id": row.matched_id,
        }


def collate_batch(items: List[Dict], descriptor_dim: int):
    batch_size = len(items)
    q_max = max(item["query_desc"].shape[0] for item in items)
    m_max = max(item["matched_desc"].shape[0] for item in items)

    q_desc = np.zeros((batch_size, q_max, descriptor_dim), dtype=np.float32)
    q_xyz = np.zeros((batch_size, q_max, 3), dtype=np.float32)
    q_mask = np.zeros((batch_size, q_max), dtype=np.float32)
    m_desc = np.zeros((batch_size, m_max, descriptor_dim), dtype=np.float32)
    m_xyz = np.zeros((batch_size, m_max, 3), dtype=np.float32)
    m_mask = np.zeros((batch_size, m_max), dtype=np.float32)
    labels = np.zeros((batch_size, 1), dtype=np.float32)
    pair_ids: List[Tuple[str, str]] = []

    for i, item in enumerate(items):
        qn = item["query_desc"].shape[0]
        mn = item["matched_desc"].shape[0]

        q_desc[i, :qn] = item["query_desc"]
        q_xyz[i, :qn] = item["query_xyz"]
        q_mask[i, :qn] = 1.0

        m_desc[i, :mn] = item["matched_desc"]
        m_xyz[i, :mn] = item["matched_xyz"]
        m_mask[i, :mn] = 1.0

        labels[i, 0] = float(item["label"])
        pair_ids.append((item["query_id"], item["matched_id"]))

    return {
        "query_desc": q_desc,
        "query_xyz": q_xyz,
        "query_mask": q_mask,
        "matched_desc": m_desc,
        "matched_xyz": m_xyz,
        "matched_mask": m_mask,
        "labels": labels,
        "pair_ids": pair_ids,
    }


def batch_iterator(dataset: ProteinPairDataset, batch_size: int, descriptor_dim: int, shuffle: bool):
    indices = np.arange(len(dataset))
    if shuffle:
        dataset.rng.shuffle(indices)
    for start in range(0, len(indices), batch_size):
        batch_idx = indices[start : start + batch_size]
        items = [dataset[int(i)] for i in batch_idx]
        yield collate_batch(items, descriptor_dim=descriptor_dim)


def sanity_check_batch(batch: Dict):
    for key in ("query_desc", "query_xyz", "matched_desc", "matched_xyz", "labels"):
        arr = batch[key]
        if not np.isfinite(arr).all():
            raise ValueError(f"Non-finite values detected in {key}")
    if batch["query_desc"].shape[0] != batch["labels"].shape[0]:
        raise ValueError("Batch size mismatch between features and labels.")

