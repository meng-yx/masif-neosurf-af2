#!/usr/bin/env python3
"""
preprocess_manifest.py — run MaSIF preprocessing for each row in an input CSV subset.

Reads pdb_path, target, ligand, ligand_path; preserves all input columns and appends
status (success | skipped | error) and error_message.

Usage:
    python scripts/python/preprocess_manifest.py \
        --input_csv data/processing/2_masif_preprocess/input_subsets/input_1.csv \
        --out_csv   data/processing/2_masif_preprocess/output_subsets/output_1.csv \
        [--preprocess_root data/preprocess] \
        [--overwrite]
"""

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from tqdm import tqdm

REQUIRED_COLUMNS = {"pdb_path", "target", "ligand", "ligand_path"}
DESC_SUBDIR = "descriptors/sc05/all_feat"
STATUS_COL = "status"
ERROR_COL = "error_message"


def descriptor_outputs_exist(output_root, target):
    out_desc_dir = Path(output_root) / DESC_SUBDIR / str(target)
    return (
        (out_desc_dir / "p1_desc_straight.npy").is_file()
        and (out_desc_dir / "p1_desc_flipped.npy").is_file()
    )


def _extract_error_message(stderr, stdout):
    combined = f"{stderr or ''}\n{stdout or ''}"
    lines = [line.strip() for line in combined.splitlines() if line.strip()]
    error_prefixes = (
        "AssertionError",
        "ValueError",
        "RuntimeError",
        "KeyError",
        "FileNotFoundError",
        "Boost.Python.ArgumentError",
    )
    for line in reversed(lines):
        if line.startswith(error_prefixes):
            return line[:500]
        match = re.search(
            r"\b(?:Assertion|Value|Runtime|Key|FileNotFound|Type|Index)Error:.*",
            line,
        )
        if match:
            return match.group(0)[:500]
        if "ArgumentError" in line:
            return line[:500]
    for line in reversed(lines):
        if "error" in line.lower() or "Traceback" in line:
            return line[:500]
    return (lines[-1][:500] if lines else "non-zero exit code")


def run_preprocess_row(
    pdb_path, target, ligand, ligand_path, output_root, repo_root, overwrite=False
):
    if not overwrite and descriptor_outputs_exist(output_root, target):
        return "skipped", ""

    env = os.environ.copy()
    if overwrite:
        env["PREPROCESS_OVERWRITE"] = "1"
    else:
        env.pop("PREPROCESS_OVERWRITE", None)

    cmd = [
        "bash",
        "scripts/bash/preprocess.sh",
        str(pdb_path),
        str(target),
        str(ligand),
        str(ligand_path),
        str(output_root),
    ]
    result = subprocess.run(
        cmd,
        cwd=repo_root,
        env=env,
        capture_output=True,
        text=True,
    )
    if result.returncode == 0:
        return "success", ""
    return "error", _extract_error_message(result.stderr, result.stdout)


def main():
    parser = argparse.ArgumentParser(
        description="Run MaSIF preprocessing for each row in a CSV subset.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input_csv", required=True, help="Input subset CSV.")
    parser.add_argument("--out_csv", required=True, help="Output CSV with status columns.")
    parser.add_argument(
        "--preprocess_root",
        default="data/preprocess",
        help="MaSIF preprocess output root (TARGET_PREPROCESS_ROOT).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        default=False,
        help="Re-run preprocessing even when descriptor files already exist.",
    )
    parser.add_argument(
        "--repo_root",
        default=".",
        help="Repository root (cwd for preprocess.sh).",
    )
    args = parser.parse_args()

    repo_root = Path(args.repo_root).resolve()
    if not (repo_root / "scripts/bash/preprocess.sh").is_file():
        print(
            f"Error: preprocess.sh not found under repo root {repo_root}",
            file=sys.stderr,
        )
        sys.exit(1)

    df = pd.read_csv(args.input_csv)
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        print(
            f"Error: --input_csv is missing required columns: {sorted(missing)}",
            file=sys.stderr,
        )
        sys.exit(1)

    df[STATUS_COL] = None
    df[ERROR_COL] = None

    n_failed = 0
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Preprocessing"):
        target = str(row.target)
        try:
            status, error_message = run_preprocess_row(
                pdb_path=row.pdb_path,
                target=target,
                ligand=row.ligand,
                ligand_path=row.ligand_path,
                output_root=args.preprocess_root,
                repo_root=repo_root,
                overwrite=args.overwrite,
            )
        except Exception as exc:
            status, error_message = "error", str(exc)[:500]

        df.at[idx, STATUS_COL] = status
        df.at[idx, ERROR_COL] = error_message
        if status == "error":
            n_failed += 1
            tqdm.write(f"Error processing {target}: {error_message}", file=sys.stderr)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)

    n_ok = len(df) - n_failed
    print(f"Wrote {out_csv}  ({n_ok}/{len(df)} rows succeeded or skipped)")
    if n_failed:
        print(f"{n_failed} row(s) failed.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
