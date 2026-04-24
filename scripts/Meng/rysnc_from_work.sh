#!/bin/bash
#SBATCH --job-name=rsync
#SBATCH --partition=standard
#SBATCH --output=/scratch/%u/logs/rsync-%A.out  # Fixed typo in filename
#SBATCH --error=/scratch/%u/logs/rsync-%A.out   # Fixed typo in filename
#SBATCH --time=08:00:00  # Reduced time limit for rsync
#SBATCH --cpus-per-task=8
#SBATCH --mem=56000

set -euo pipefail

if ! command -v pigz >/dev/null 2>&1; then
  echo "Error: pigz is required but not found in PATH." >&2
  exit 1
fi


# Source directory to be copied from
DEST_DIR=$(git rev-parse --show-toplevel)
DEST_DIR_NAME=$(basename "$DEST_DIR")
DEST_PARENT_DIR=$(dirname "$DEST_DIR")
TAR_NAME="${DEST_DIR_NAME}.tar.gz"
LOCAL_TAR_PATH="${DEST_PARENT_DIR}/${TAR_NAME}"
EXTRACT_ROOT="${DEST_PARENT_DIR}/${DEST_DIR_NAME}_from_work_extract"
EXTRACTED_DIR="${EXTRACT_ROOT}/${DEST_DIR_NAME}"

# Destination to be copied towards
SRC_TAR_PATH="/work/upthomae/Meng/JED_TO_KUMA/${TAR_NAME}"

echo "SRC_TAR_PATH: $SRC_TAR_PATH"
echo "DEST_DIR: $DEST_DIR"
echo "Starting tarball pull at $(date)"
rsync -avh --progress "$SRC_TAR_PATH" "$LOCAL_TAR_PATH"

echo "Extracting tarball at $(date)"
rm -rf "$EXTRACT_ROOT"
mkdir -p "$EXTRACT_ROOT"
tar --use-compress-program=pigz -xf "$LOCAL_TAR_PATH" -C "$EXTRACT_ROOT"

echo "Syncing extracted directory into DEST_DIR at $(date)"
rsync -avh --progress "$EXTRACTED_DIR/" "$DEST_DIR/"
echo "Sync from tarball completed at $(date)"
