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
SRC_DIR=$(git rev-parse --show-toplevel)
SRC_DIR_NAME=$(basename "$SRC_DIR")
SRC_PARENT_DIR=$(dirname "$SRC_DIR")
TAR_NAME="${SRC_DIR_NAME}.tar.gz"
TAR_PATH="${SRC_PARENT_DIR}/${TAR_NAME}"

# Destination to be copied towards
DEST_DIR="/work/upthomae/Meng/JED_TO_KUMA"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

echo "Creating tarball at $(date)"
tar -C "$SRC_PARENT_DIR" --use-compress-program=pigz -cf "$TAR_PATH" "$SRC_DIR_NAME"

echo "Starting rsync at $(date)"
rsync -avh --progress "$TAR_PATH" "$DEST_DIR/"
echo "Tarball rsync completed at $(date)"