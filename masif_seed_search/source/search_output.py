"""
CSV output for MaSIF search hits (one file per matched protein per site).
"""
import csv
import os
from pathlib import Path

import numpy as np
import pandas as pd

HIT_CSV_COLUMNS = [
    "target",
    "target_site",
    "target_vix",
    "matched_protein",
    "matched_patch_id",
    "score",
    "desc_dist_score",
    "clashing_ca",
    "clashing_heavy",
    "matched_vix",
    "desc_dist",
    "iface_score",
    "mean_desc_dist_score",
    "flattened_transform",
]

CLUSTER_COLUMNS = ["cluster_id", "cluster_size", "cluster_mean_rmsd"]
CLUSTERED_HIT_CSV_COLUMNS = HIT_CSV_COLUMNS + CLUSTER_COLUMNS


def hit_csv_path(site_outdir, matched_protein):
    return os.path.join(site_outdir, f"{matched_protein}.csv")


# ---------------------------------------------------------------------------
# Site-level completion flags  (single source of truth for resume)
# ---------------------------------------------------------------------------

def site_completed_flag_path(site_outdir, progress_id):
    """Return path to the per-task site-completion sentinel file."""
    return os.path.join(site_outdir, ".progress", f"SUBSET_{progress_id}_COMPLETED")


def is_site_completed(site_outdir, progress_id):
    """Return True when the SUBSET_<progress_id>_COMPLETED flag exists for this site."""
    if progress_id is None:
        return False
    return os.path.isfile(site_completed_flag_path(site_outdir, progress_id))


def mark_site_completed(site_outdir, progress_id):
    """Write the SUBSET_<progress_id>_COMPLETED sentinel file for this site."""
    if progress_id is None:
        return
    prog_dir = os.path.join(site_outdir, ".progress")
    os.makedirs(prog_dir, exist_ok=True)
    flag = site_completed_flag_path(site_outdir, progress_id)
    with open(flag, "w") as f:
        f.write("")
        f.flush()
        os.fsync(f.fileno())


def build_hit_row(
    params,
    matched_protein,
    matched_patch_id,
    score,
    desc_dist_score,
    mean_desc_dist_score,
    clashing_ca,
    clashing_heavy,
    matched_vix,
    transformation,
    first_stage_scores=None,
    patch_index=0,
):
    """Build one CSV row dict for a passing alignment hit."""
    row = {
        "target": params["target_name"],
        "target_site": params["target_site"],
        "target_vix": int(params["target_vix"]),
        "matched_protein": matched_protein,
        "matched_patch_id": int(matched_patch_id),
        "score": float(score),
        "desc_dist_score": float(desc_dist_score),
        "clashing_ca": int(clashing_ca),
        "clashing_heavy": int(clashing_heavy),
        "matched_vix": int(matched_vix),
        "mean_desc_dist_score": float(mean_desc_dist_score),
        "flattened_transform": ",".join(map(str, np.asarray(transformation).reshape(4, 4).flatten())),
    }
    if first_stage_scores is not None:
        row["desc_dist"] = float(first_stage_scores["desc_dist"][patch_index])
        row["iface_score"] = float(first_stage_scores["iface_score"][patch_index])
    else:
        row["desc_dist"] = ""
        row["iface_score"] = ""
    return row


def write_site_seed_hits(site_outdir, matched_protein, rows, progress_id=None):
    """Write hit CSV for one (site, seed). Does nothing when rows is empty."""
    if not rows:
        return
    os.makedirs(site_outdir, exist_ok=True)
    csv_path = hit_csv_path(site_outdir, matched_protein)
    tmp_path = f"{csv_path}.tmp"
    with open(tmp_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=HIT_CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in HIT_CSV_COLUMNS})
    os.replace(tmp_path, csv_path)


def gather_search_results(results_root, query_target):
    """Concatenate per-protein hit CSVs under results_root/query_target/site_*."""
    target_dir = Path(results_root) / query_target
    if not target_dir.is_dir():
        raise FileNotFoundError(f"No search output directory: {target_dir}")

    frames = []
    for site_dir in sorted(target_dir.glob("site_*")):
        for csv_path in sorted(site_dir.glob("*.csv")):
            part = pd.read_csv(csv_path)
            if len(part) > 0:
                frames.append(part)

    if not frames:
        return pd.DataFrame(columns=HIT_CSV_COLUMNS)
    return pd.concat(frames, ignore_index=True)


def gather_clustered_results(results_root, query_target):
    """Concatenate per-seed clustered CSVs under results_root/query_target/clustered_matches/."""
    clustered_dir = Path(results_root) / query_target / "clustered_matches"
    if not clustered_dir.is_dir():
        return pd.DataFrame(columns=CLUSTERED_HIT_CSV_COLUMNS)

    frames = []
    for csv_path in sorted(clustered_dir.glob("*.csv")):
        part = pd.read_csv(csv_path)
        if len(part) > 0:
            frames.append(part)

    if not frames:
        return pd.DataFrame(columns=CLUSTERED_HIT_CSV_COLUMNS)
    return pd.concat(frames, ignore_index=True)
