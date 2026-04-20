#!/bin/bash
#SBATCH --job-name=rsync
#SBATCH --partition=standard
#SBATCH --output=logs/rsync-%A.out  # Fixed typo in filename
#SBATCH --error=logs/rsync-%A.out   # Fixed typo in filename
#SBATCH --time=02:00:00  # Reduced time limit for rsync
#SBATCH --mem=7000


# Source directory to be copied from
DEST_DIR=$(git rev-parse --show-toplevel)

DEST_DIR_NAME=$(basename "$DEST_DIR")

# Destination to be copied towards
SRC_DIR="/work/upthomae/Meng/JED_TO_KUMA/$DEST_DIR_NAME/"

echo "SRC_DIR: $SRC_DIR"
echo "DEST_DIR: $DEST_DIR"
echo "Starting rsync at $(date)"
rsync -avh --progress "$SRC_DIR" "$DEST_DIR"
echo "Rsync completed at $(date)"


