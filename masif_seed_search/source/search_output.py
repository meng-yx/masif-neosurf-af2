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
# Progress-manifest helpers (zero-hit resume without header-only CSVs)
# ---------------------------------------------------------------------------

def progress_file_path(site_outdir, progress_id):
    """Return the path to the per-task progress file inside a site directory."""
    return os.path.join(site_outdir, ".progress", f"{progress_id}.txt")


def load_completed_seeds(site_outdir, progress_id):
    """
    Return the set of seed IDs recorded as complete in the progress manifest.
    Missing file is treated as empty (no seeds done yet).
    """
    p = progress_file_path(site_outdir, progress_id)
    if not os.path.isfile(p):
        return set()
    with open(p) as f:
        return {line.strip() for line in f if line.strip()}


def mark_site_seed_complete(site_outdir, progress_id, seed_id):
    """
    Append seed_id to the per-task progress file.  Uses flush+fsync so the
    entry survives an abrupt SLURM job timeout.
    """
    prog_dir = os.path.join(site_outdir, ".progress")
    os.makedirs(prog_dir, exist_ok=True)
    with open(os.path.join(prog_dir, f"{progress_id}.txt"), "a") as f:
        f.write(f"{seed_id}\n")
        f.flush()
        os.fsync(f.fileno())


def filter_seeds_for_resume(site_outdir, seed_ids, resume, progress_id=None):
    """
    Return the subset of seed_ids that still need to be processed.

    In resume mode a seed is considered complete when:
    - its hit CSV exists on disk (new hits, or legacy header-only CSV), OR
    - its ID appears in the per-task progress manifest (zero-hit, new scheme).
    """
    if not resume:
        return list(seed_ids)

    # Load progress manifest once for all seeds at this site.
    completed = load_completed_seeds(site_outdir, progress_id) if progress_id else set()

    pending = []
    for seed_id in seed_ids:
        if seed_id in completed:
            print(f"Resume: skipping {seed_id} at {os.path.basename(site_outdir)} "
                  f"(found in progress manifest)")
            continue
        if os.path.isfile(hit_csv_path(site_outdir, seed_id)):
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


def write_site_seed_hits(site_outdir, matched_protein, rows, progress_id=None):
    """
    Write all hits for one (site, seed) in a single atomic CSV update.

    When rows is non-empty the hit CSV is written as before.  When rows is
    empty, no CSV is written; instead the seed is recorded in the per-task
    progress manifest (requires progress_id).  Legacy callers that omit
    progress_id still get a header-only CSV for backward compatibility.
    """
    os.makedirs(site_outdir, exist_ok=True)
    if rows:
        csv_path = hit_csv_path(site_outdir, matched_protein)
        tmp_path = f"{csv_path}.tmp"
        with open(tmp_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=HIT_CSV_COLUMNS, extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow({col: row.get(col, "") for col in HIT_CSV_COLUMNS})
        os.replace(tmp_path, csv_path)
    elif progress_id is not None:
        mark_site_seed_complete(site_outdir, progress_id, matched_protein)
    else:
        # Legacy fallback: header-only CSV keeps old resume logic working.
        csv_path = hit_csv_path(site_outdir, matched_protein)
        tmp_path = f"{csv_path}.tmp"
        with open(tmp_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=HIT_CSV_COLUMNS, extrasaction="ignore")
            writer.writeheader()
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
