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


def hit_csv_path(site_outdir, matched_protein):
    return os.path.join(site_outdir, f"{matched_protein}.csv")


def is_site_seed_complete(site_outdir, matched_protein):
    """A CSV on disk means align_protein finished for this (site, seed)."""
    return os.path.isfile(hit_csv_path(site_outdir, matched_protein))


def filter_seeds_for_resume(site_outdir, seed_ids, resume):
    if not resume:
        return list(seed_ids)
    pending = []
    for seed_id in seed_ids:
        if is_site_seed_complete(site_outdir, seed_id):
            print(f"Resume: skipping {seed_id} at {os.path.basename(site_outdir)} "
                  f"(found {hit_csv_path(site_outdir, seed_id)})")
            continue
        pending.append(seed_id)
    return pending


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


def write_site_seed_hits(site_outdir, matched_protein, rows):
    """
    Write all hits for one (site, seed) in a single atomic CSV update.
    Header-only file means the search completed with zero passing hits.
    """
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
            frames.append(pd.read_csv(csv_path))

    if not frames:
        return pd.DataFrame(columns=HIT_CSV_COLUMNS)
    return pd.concat(frames, ignore_index=True)


def gather_clustered_results(results_root, query_target):
    """Concatenate per-seed clustered CSVs under results_root/query_target/clustered_matches/."""
    clustered_dir = Path(results_root) / query_target / "clustered_matches"
    if not clustered_dir.is_dir():
        raise FileNotFoundError(f"No clustered output directory: {clustered_dir}")

    frames = []
    for csv_path in sorted(clustered_dir.glob("*.csv")):
        frames.append(pd.read_csv(csv_path))

    if not frames:
        return pd.DataFrame(columns=HIT_CSV_COLUMNS)
    return pd.concat(frames, ignore_index=True)
