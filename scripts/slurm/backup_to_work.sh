#!/bin/bash
#SBATCH --job-name=backup_to_work
#SBATCH --output=logs/backup_to_work-%j.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14000

: '
Back up the updated (gitignored) result data from this scratch clone to the /work
backup. Mirror image of restore_from_work.sh.

Submit from repo root:

  sbatch scripts/slurm/backup_to_work.sh                    # processing/ + masif_search tar
  sbatch scripts/slurm/backup_to_work.sh --dry-run          # preview, no writes
  sbatch scripts/slurm/backup_to_work.sh --with-preprocess  # also rsync the ~105 GB preprocess tree

Backs up:
  data/processing   -> /work/.../processing                                    (result tables)
  data/masif_search -> /work/.../human_reference_.../masif_search.tar.gz       (re-tarred, atomic)
  data/preprocess   -> /work/.../preprocess          (only with --with-preprocess)

restore_from_work.sh reads masif_search back from that tar.gz, so this keeps the
round trip working. rsync is append/update only (no --delete): stale files on the
backup are left in place.
'

set -euo pipefail

DEST_ROOT=/work/upthomae/Meng/Neosurf_Neosurf
SRC_ROOT="$(pwd)/data"
DRY_RUN=0
WITH_PREPROCESS=0

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        --with-preprocess) WITH_PREPROCESS=1 ;;
        -h|--help) sed -n "2,25p" "$0"; exit 0 ;;
        *) echo "Unknown argument: $arg" >&2
           echo "Usage: sbatch $0 [--dry-run] [--with-preprocess]" >&2; exit 1 ;;
    esac
done

if [[ ! -f scripts/config.sh ]]; then
    echo "Error: run sbatch from repo root (scripts/config.sh not found in $(pwd))" >&2
    exit 1
fi

mkdir -p logs
[[ -d "$DEST_ROOT" ]] || { echo "Error: work backup root missing: $DEST_ROOT" >&2; exit 1; }

RSYNC_FLAGS=(-aH --info=progress2 --partial)
[[ "$DRY_RUN" -eq 1 ]] && RSYNC_FLAGS+=(--dry-run)

echo "Job ID:     ${SLURM_JOB_ID:-local}"
echo "Host:       $(hostname)"
echo "Source:     $SRC_ROOT"
echo "Dest:       $DEST_ROOT"
echo "Dry run:    $DRY_RUN"
echo "Preprocess: $WITH_PREPROCESS"
echo "Started:    $(date)"
echo
df -h "$DEST_ROOT" "$SRC_ROOT" 2>/dev/null || true
echo

fail=0

# 1) processing/ -> work (result tables: 3_masif_search + 4_enrich_metrics)
echo "================================================================"
echo "rsync  $SRC_ROOT/processing/  ->  $DEST_ROOT/processing/"
echo "size:  $(du -sh "$SRC_ROOT/processing" 2>/dev/null | awk '{print $1}')"
echo "================================================================"
mkdir -p "$DEST_ROOT/processing"
if rsync "${RSYNC_FLAGS[@]}" "$SRC_ROOT/processing/" "$DEST_ROOT/processing/"; then
    echo "OK: processing"
else
    echo "FAILED: rsync processing (exit $?)" >&2; fail=1
fi
echo

# 2) masif_search -> re-tar to human_reference_.../masif_search.tar.gz (atomic)
TAR_DEST="$DEST_ROOT/human_reference_proteome_liganded_pdbs_260529/masif_search.tar.gz"
echo "================================================================"
echo "re-tar  $SRC_ROOT/masif_search  ->  $TAR_DEST"
echo "================================================================"
if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[dry-run] would tar data/masif_search -> $TAR_DEST"
else
    mkdir -p "$(dirname "$TAR_DEST")"
    tmp="${TAR_DEST}.tmp.${SLURM_JOB_ID:-$$}"
    # Entries are prefixed masif_search/ so restore_from_work.sh extracts correctly.
    if tar -czf "$tmp" -C "$SRC_ROOT" masif_search; then
        mv -f "$tmp" "$TAR_DEST"
        echo "OK: wrote $TAR_DEST ($(du -sh "$TAR_DEST" 2>/dev/null | awk '{print $1}'))"
    else
        echo "FAILED: tar masif_search" >&2; rm -f "$tmp"; fail=1
    fi
fi
echo

# 3) preprocess (~105 GB) -- only when explicitly requested (unchanged by the search stage)
if [[ "$WITH_PREPROCESS" -eq 1 ]]; then
    echo "================================================================"
    echo "rsync  $SRC_ROOT/preprocess/  ->  $DEST_ROOT/preprocess/   (large)"
    echo "================================================================"
    mkdir -p "$DEST_ROOT/preprocess"
    if rsync "${RSYNC_FLAGS[@]}" "$SRC_ROOT/preprocess/" "$DEST_ROOT/preprocess/"; then
        echo "OK: preprocess"
    else
        echo "FAILED: rsync preprocess (exit $?)" >&2; fail=1
    fi
    echo
fi

echo "================================================================"
echo "Backup summary (work)"
echo "================================================================"
du -sh "$DEST_ROOT"/processing "$TAR_DEST" 2>/dev/null || true
echo
echo "Finished:  $(date)"
echo "Exit code: $fail"
exit "$fail"
