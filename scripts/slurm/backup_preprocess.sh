#!/bin/bash
#SBATCH --job-name=backup_preprocess
#SBATCH --output=logs/backup_preprocess-%j.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=14000

: '
Back up search-required files from a MaSIF preprocess tree.

Submit from repo root:

  # Copy (~100 GB; use before deleting the full preprocess tree)
  sbatch scripts/slurm/backup_preprocess.sh rsync data/preprocess data/preprocess_search

  # Symlink mirror (fast; original must remain for links to resolve)
  sbatch scripts/slurm/backup_preprocess.sh symlink data/preprocess data/preprocess_search

  # Dry run
  sbatch scripts/slurm/backup_preprocess.sh rsync data/preprocess /tmp/preprocess_search --dry-run

Arguments:
  $1  MODE   rsync | symlink  (default: rsync)
  $2  SOURCE preprocess root (default: data/preprocess)
  $3  DEST   output directory (required)
  $4  optional: --dry-run
'

MODE=symlink
SOURCE=data/preprocess
DEST=data/preprocess_backup
EXTRA=""


mkdir -p logs

echo "Job ID:    ${SLURM_JOB_ID:-local}"
echo "Mode:      ${MODE}"
echo "Source:    ${SOURCE}"
echo "Dest:      ${DEST}"
echo "Extra:     ${EXTRA}"
echo "Started:   $(date)"
echo

BACKUP_ARGS=("${MODE}" "${SOURCE}" "${DEST}")
if [[ -n "$EXTRA" ]]; then
    BACKUP_ARGS+=("${EXTRA}")
fi

bash scripts/bash/backup_masif_preprocess_dir.sh "${BACKUP_ARGS[@]}"
status=$?

echo
echo "Finished:  $(date)"
echo "Exit code: ${status}"
exit "${status}"
