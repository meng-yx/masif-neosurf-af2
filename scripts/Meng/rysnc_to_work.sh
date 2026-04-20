#!/bin/bash
#SBATCH --job-name=rsync
#SBATCH --partition=standard
#SBATCH --output=logs/rsync-%A.out  # Fixed typo in filename
#SBATCH --error=logs/rsync-%A.out   # Fixed typo in filename
#SBATCH --time=02:00:00  # Reduced time limit for rsync
#SBATCH --mem=7000


# Source directory to be copied from
SRC_DIR=$(git rev-parse --show-toplevel)

SRC_DIR_NAME=$(basename "$SRC_DIR")

# Destination to be copied towards
DEST_DIR="/work/upthomae/Meng/JED_TO_KUMA/$SRC_DIR_NAME"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

echo "Starting rsync at $(date)"
rsync -avh --progress "$SRC_DIR/" "$DEST_DIR/"
echo "Rsync completed at $(date)"