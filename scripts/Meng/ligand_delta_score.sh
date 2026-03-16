#!/bin/bash
#SBATCH --job-name=ligand-delta
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ligand-delta-%A/slurm-%A_%a.out
#SBATCH --error=logs/ligand-delta-%A/slurm-%A_%a.out

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: sbatch --array 1-<n_jobs> ligand_delta_score.sh <input_dir> <output_dir>" >&2
  exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID is not set. Run this script with sbatch --array." >&2
  exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
TASK_ID="${SLURM_ARRAY_TASK_ID}"
INPUT_CSV="${INPUT_DIR}/input_${TASK_ID}.csv"
OUTPUT_CSV="${OUTPUT_DIR}/output_${TASK_ID}.csv"
FAILED_DIR="${OUTPUT_DIR}/failed"
mkdir -p "${FAILED_DIR}"
FAILED_CSV="${FAILED_DIR}/output_${TASK_ID}.failed.csv"

# Configurable runtime parameters.
MASIF_NEOSURF_REPO="${MASIF_NEOSURF_REPO:-/scratch/ymeng/masif-neosurf-af2}"
PYTHON_BIN="${PYTHON_BIN:-python}"
IMAGE="${IMAGE:-/scratch/ymeng/tricomplex-design_Meng/masif-neosurf_v0.1.sif}"
SINGULARITY_BIND="${SINGULARITY_BIND:-${MASIF_NEOSURF_REPO}:${MASIF_NEOSURF_REPO},/work:/work,/scratch:/scratch}"
: "${TARGET_PREPROC_DIR:?Set TARGET_PREPROC_DIR (e.g. /path/to/target/preprocessing)}"
: "${MATCH_PREPROC_DIR:?Set MATCH_PREPROC_DIR (e.g. /path/to/match/preprocessing)}"

mkdir -p "${OUTPUT_DIR}" logs

if [[ ! -f "${INPUT_CSV}" ]]; then
  echo "ERROR: Missing input file: ${INPUT_CSV}" >&2
  exit 1
fi

echo "SLURM_JOB_ID=${SLURM_JOB_ID:-unknown} TASK_ID=${TASK_ID}"
echo "Input:  ${INPUT_CSV}"
echo "Output: ${OUTPUT_CSV}"
echo "Repo:   ${MASIF_NEOSURF_REPO}"
echo "Image:  ${IMAGE}"

# Ensure required Python deps exist inside the container runtime.
singularity exec --cleanenv --bind "${SINGULARITY_BIND}" "${IMAGE}" \
  bash -lc "pip install pandas numpy tqdm"

singularity exec --cleanenv --bind "${SINGULARITY_BIND}" "${IMAGE}" \
  "${PYTHON_BIN}" - "${INPUT_CSV}" "${OUTPUT_CSV}" "${FAILED_CSV}" "${MASIF_NEOSURF_REPO}" "${TARGET_PREPROC_DIR}" "${MATCH_PREPROC_DIR}" "${PYTHON_BIN}" "${TASK_ID}" <<'PY'
import csv
import json
import pathlib
import subprocess
import sys
from typing import Dict, List
from tqdm import tqdm

input_csv = pathlib.Path(sys.argv[1])
output_csv = pathlib.Path(sys.argv[2])
failed_csv = pathlib.Path(sys.argv[3])
repo = pathlib.Path(sys.argv[4])
target_preproc_dir = pathlib.Path(sys.argv[5])
match_preproc_dir = pathlib.Path(sys.argv[6])
python_bin = sys.argv[7]
task_id = sys.argv[8]


def parse_json_blob(text: str) -> Dict[str, object]:
    start = text.find("{")
    end = text.rfind("}")
    if start == -1 or end == -1 or end < start:
        raise ValueError("No JSON object found in command output.")
    return json.loads(text[start : end + 1])


def run_one(row: Dict[str, str]) -> Dict[str, object]:
    cmd = [
        python_bin,
        str(repo / "ligand_delta_scores.py"),
        "--target_preproc_dir",
        str(target_preproc_dir),
        "--target_name",
        str(row["target"]),
        "--match_preproc_dir",
        str(match_preproc_dir),
        "--match_name",
        str(row["matched_protein"]),
        "--match_vix",
        str(int(float(row["matched_vix"]))),
        f"--flattened_transform={row['flattened_transform']}",
        "--masif-neosurf-repo",
        str(repo),
    ]

    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or proc.stdout.strip() or "Unknown error.")
    return parse_json_blob(proc.stdout)


with input_csv.open("r", newline="") as handle:
    reader = csv.DictReader(handle)
    input_rows: List[Dict[str, str]] = list(reader)

if not input_rows:
    # Keep behavior explicit for empty array chunks.
    with output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["status", "message"])
        writer.writerow(["empty", "No rows in this split file"])
    print(f"Input file is empty: {input_csv}")
    sys.exit(0)

ok_rows: List[Dict[str, object]] = []
failed_rows: List[Dict[str, object]] = []

rows_iter = tqdm(
    input_rows,
    total=len(input_rows),
    desc=f"Task {task_id}",
    unit="row",
)
for idx, row in enumerate(rows_iter, start=1):
    try:
        result = run_one(row)
        merged = dict(row)
        merged.update(result)
        ok_rows.append(merged)
    except Exception as exc:  # pylint: disable=broad-except
        fail_row = dict(row)
        fail_row["error"] = str(exc)
        fail_row["row_index_1based"] = idx
        failed_rows.append(fail_row)

if ok_rows:
    fieldnames: List[str] = []
    for r in ok_rows:
        for key in r.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(ok_rows)

if failed_rows:
    fail_fields: List[str] = []
    for r in failed_rows:
        for key in r.keys():
            if key not in fail_fields:
                fail_fields.append(key)
    with failed_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fail_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(failed_rows)

print(f"Processed rows: {len(input_rows)}")
print(f"Succeeded: {len(ok_rows)} -> {output_csv}")
print(f"Failed: {len(failed_rows)}" + (f" -> {failed_csv}" if failed_rows else ""))

if failed_rows:
    sys.exit(2)
PY
