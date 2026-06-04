"""
Structural clustering of MaSIF search hits (DBSCAN on backbone RMSD).
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from sklearn.cluster import DBSCAN

from search_output import HIT_CSV_COLUMNS

RMSD_THRESH = 5.0
CLUSTER_COLUMNS = ["cluster_id", "cluster_size", "cluster_mean_rmsd"]


def get_transformed_struct(row, pdb_dir):
    """Apply MaSIF alignment transform to the preprocessed seed structure."""
    pdb_path = Path(pdb_dir) / f"{row.matched_protein}.pdb"
    struct = PDBParser(QUIET=True).get_structure(row.matched_protein, pdb_path)
    transform = np.array(list(map(float, row.flattened_transform.split(",")))).reshape(4, 4)
    struct.transform(rot=transform[:3, :3].T, tran=transform[:3, 3])
    return struct


def structural_clusters(df, rmsd_thresh=RMSD_THRESH, pdb_dir=None):
    """
    Cluster alignments per matched_protein by backbone RMSD (DBSCAN).
    Adds cluster_id, cluster_size, cluster_mean_rmsd.
    """
    if pdb_dir is None:
        raise ValueError("pdb_dir is required")
    if len(df) == 0:
        return df

    df = df.copy()

    def get_coords(struct, atoms_to_keep=("N", "CA", "C")):
        return np.stack([a.get_coord() for a in struct.get_atoms() if a.name in atoms_to_keep])

    def rmsd(coords1, coords2):
        assert coords1.shape == coords2.shape
        diff = coords1 - coords2
        return np.sqrt(np.mean(np.sum(diff * diff, axis=-1)))

    current_cluster_label = 0
    for domain_name in df.matched_protein.unique():
        domain_table = df[df.matched_protein == domain_name]

        if len(domain_table) == 1:
            df.loc[domain_table.index, "cluster_id"] = current_cluster_label
            df.loc[domain_table.index, "cluster_size"] = 1
            df.loc[domain_table.index, "cluster_mean_rmsd"] = 0.0
            current_cluster_label += 1
            continue

        positions = [
            get_coords(get_transformed_struct(row, pdb_dir=pdb_dir))
            for _, row in domain_table.iterrows()
        ]

        rmsd_vals = np.zeros((len(domain_table), len(domain_table)))
        for i in range(len(domain_table)):
            for j in range(i + 1, len(domain_table)):
                rmsd_ij = rmsd(positions[i], positions[j])
                rmsd_vals[i, j] = rmsd_ij
                rmsd_vals[j, i] = rmsd_ij

        labels = DBSCAN(eps=rmsd_thresh, min_samples=2, metric="precomputed").fit_predict(rmsd_vals)
        labels[labels == -1] = np.arange((labels == -1).sum()) + labels.max() + 1

        for lb in set(labels):
            mask = labels == lb
            cluster_table = domain_table.iloc[mask]
            df.loc[cluster_table.index, "cluster_id"] = current_cluster_label
            df.loc[cluster_table.index, "cluster_size"] = len(cluster_table)
            df.loc[cluster_table.index, "cluster_mean_rmsd"] = rmsd_vals[mask][:, mask].mean(axis=1)
            current_cluster_label += 1

    return df


def deduplicate_by_cluster(df, group_cols=("target", "matched_protein", "cluster_id"), score_col="score"):
    """Keep the highest-scoring row per group."""
    idx = df.groupby(list(group_cols))[score_col].idxmax()
    return df.loc[idx].reset_index(drop=True)


def gather_seed_hits_across_sites(target_dir, seed_id):
    """Concatenate hit rows for one seed from all site_* CSVs (skips missing / header-only)."""
    frames = []
    for site_dir in sorted(Path(target_dir).glob("site_*")):
        csv_path = site_dir / f"{seed_id}.csv"
        if not csv_path.is_file():
            continue
        part = pd.read_csv(csv_path)
        if len(part) > 0:
            frames.append(part)
    if not frames:
        return pd.DataFrame(columns=HIT_CSV_COLUMNS)
    return pd.concat(frames, ignore_index=True)


def clustered_csv_path(target_dir, seed_id):
    return Path(target_dir) / "clustered_matches" / f"{seed_id}.csv"


def write_clustered_hits(target_dir, seed_id, df):
    """Atomically write clustered hits for one seed."""
    out_dir = Path(target_dir) / "clustered_matches"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = clustered_csv_path(target_dir, seed_id)
    tmp_path = f"{csv_path}.tmp"
    df.to_csv(tmp_path, index=False)
    os.replace(tmp_path, csv_path)


def cluster_search_hits_for_subset(params, rmsd_thresh=RMSD_THRESH, verbose=True):
    """
    After search completes for this process, cluster each seed in seed_ppi_pair_ids
    from site_* CSVs and write clustered_matches/{seed_id}.csv.
    """
    target_dir = Path(params["out_dir"]) / params["target_name"]
    if not target_dir.is_dir():
        print(f"Clustering: no output directory {target_dir}, skipping.")
        return

    pdb_dir = params["seed_pdb_dir"]
    seed_ids = list(params["seed_ppi_pair_ids"])
    if verbose:
        try:
            from tqdm import tqdm
            seed_iter = tqdm(seed_ids, desc="cluster_search_hits")
        except ImportError:
            seed_iter = seed_ids
    else:
        seed_iter = seed_ids

    n_written = 0
    n_skipped = 0
    for seed_id in seed_iter:
        df = gather_seed_hits_across_sites(target_dir, seed_id)
        if len(df) == 0:
            n_skipped += 1
            continue
        df_clustered = structural_clusters(df, rmsd_thresh=rmsd_thresh, pdb_dir=pdb_dir)
        write_clustered_hits(target_dir, seed_id, df_clustered)
        n_written += 1

    print(
        f"Clustering done: wrote {n_written} file(s) under "
        f"{target_dir / 'clustered_matches'} ({n_skipped} seed(s) with no hits)."
    )
