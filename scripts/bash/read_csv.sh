#!/usr/bin/env bash
# Read one row from a CSV (header row 1) and set shell variables named after columns.
# Usage: source scripts/bash/read_csv.sh <csv_file> <row_number>
# Row 2 is the first data row. For SLURM arrays: row=$((SLURM_ARRAY_TASK_ID + 1)).

csv_file=$1
row=$2

IFS=, read -r -a headers < "$csv_file"
line=$(awk -v r="$row" 'NR==r {print; exit}' "$csv_file")

if [[ -z "$line" ]]; then
    echo "Error: no row ${row} in ${csv_file}" >&2
    return 1 2>/dev/null || exit 1
fi

IFS=, read -r -a values <<< "$line"

for i in "${!headers[@]}"; do
    header=$(echo "${headers[i]}" | tr -d '[:space:]' | tr -d '"')
    value=$(echo "${values[i]}" | sed 's/^ *//;s/ *$//' | tr -d '"')
    printf -v "$header" '%s' "$value"
done
