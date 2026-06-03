#!/bin/bash
#SBATCH --job-name=masif_search
#SBATCH --output=logs/search-%A/slurm-%A_%a.out
#SBATCH --time=01:30:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G

: '
# Submit from masif-neosurf repo root:
#   cd /path/to/masif-neosurf
#   sbatch --array=1-N scripts/slurm/search_array.sh data/masif_search/search_manifest.csv
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

for var in query_target target_chain target_resid seed_target seed_chain seed_resid out_dir; do
    if [[ -z "${!var}" ]]; then
        echo "Error: ${var} empty after reading row ${CSV_ROW} of ${CSV_FILE}" >&2
        exit 1
    fi
done

echo "Array task ${SLURM_ARRAY_TASK_ID} (CSV row ${CSV_ROW}): ${query_target} vs ${seed_target}"

bash scripts/bash/search.sh \
    "${query_target}" "${target_chain}" "${target_resid}" \
    "${seed_target}" "${seed_chain}" "${seed_resid}" "${out_dir}"
