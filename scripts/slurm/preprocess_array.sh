#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --output=logs/preprocess-%A/slurm-%A_%a.out
#SBATCH --time=0:10:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=14000

: '
# Submit from masif-neosurf repo root (Slurm job cwd = directory where you run sbatch):
#   cd /path/to/masif-neosurf
#   sbatch --array=1-28 scripts/slurm/preprocess_array.sh data/input/input_manifest.csv
'

CSV_FILE=$1
if [[ -z "$CSV_FILE" ]]; then
    echo "Usage: sbatch --array=1-N $0 <CSV_FILE>" >&2
    exit 1
fi
CSV_FILE="${CSV_FILE#\$}"
if [[ ! -f "$CSV_FILE" ]]; then
    echo "Error: CSV file not found: $CSV_FILE" >&2
    exit 1
fi

if [[ ! -f scripts/config.sh ]]; then
    echo "Error: run sbatch from masif-neosurf repo root (scripts/config.sh not found in $(pwd))" >&2
    exit 1
fi

source scripts/config.sh

: "${SLURM_ARRAY_TASK_ID:=1}"
CSV_ROW=$((SLURM_ARRAY_TASK_ID + 1))

source scripts/bash/read_csv.sh "${CSV_FILE}" "${CSV_ROW}"

PDB_PATH=$pdb_path
TARGET=$target
LIGAND=$ligand
LIGAND_PATH=$ligand_path

for var in PDB_PATH TARGET LIGAND LIGAND_PATH; do
    if [[ -z "${!var}" ]]; then
        echo "Error: ${var} empty after reading row ${CSV_ROW} of ${CSV_FILE}" >&2
        exit 1
    fi
done

echo "Array task ${SLURM_ARRAY_TASK_ID} (CSV row ${CSV_ROW}): TARGET=${TARGET}"

bash scripts/bash/preprocess.sh \
    "$PDB_PATH" "$TARGET" "$LIGAND" "$LIGAND_PATH" "$DATABASE_ROOT"
