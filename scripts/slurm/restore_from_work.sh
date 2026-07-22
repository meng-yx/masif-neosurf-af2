#!/bin/bash
#SBATCH --job-name=restore_from_work
#SBATCH --output=logs/restore_from_work-%j.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=14000

: '
Restore gitignored data from the /work backup into this scratch clone.

Submit from repo root:

  sbatch scripts/slurm/restore_from_work.sh

  # Dry run (lists what would be copied, no writes)
  sbatch scripts/slurm/restore_from_work.sh --dry-run

  # Skip extracting masif_search.tar.gz after copy
  sbatch scripts/slurm/restore_from_work.sh --no-extract

Work backup layout (top-level dirs) is mirrored under data/:

  /work/.../preprocess          -> data/preprocess          (~105 GB)
  /work/.../processing          -> data/processing          (~111 MB)
  /work/.../human_reference_... -> data/human_reference_... (~113 MB; includes masif_search.tar.gz)
  /work/.../preprocess_test     -> data/preprocess_test     (empty test dir)

After rsync, masif_search.tar.gz is extracted to data/masif_search/ unless --no-extract.
'

set -euo pipefail

SRC_ROOT=/work/upthomae/Meng/Neosurf_Neosurf
DEST_ROOT="$(pwd)/data"
DRY_RUN=0
EXTRACT_TAR=1

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=1 ;;
        --no-extract) EXTRACT_TAR=0 ;;
        -h|--help)
            sed -n "2,30p" "$0"
            exit 0
            ;;
        *)
            echo "Unknown argument: $arg" >&2
            echo "Usage: sbatch $0 [--dry-run] [--no-extract]" >&2
            exit 1
            ;;
    esac
done

if [[ ! -f scripts/config.sh ]]; then
    echo "Error: run sbatch from repo root (scripts/config.sh not found in $(pwd))" >&2
    exit 1
fi

mkdir -p logs
mkdir -p "$DEST_ROOT"

# name_on_work -> name_under_data/
RESTORE_DIRS=(
    "human_reference_proteome_liganded_pdbs_260529:human_reference_proteome_liganded_pdbs_260529"
    "preprocess:preprocess"
    "processing:processing"
    "preprocess_test:preprocess_test"
)

RSYNC_FLAGS=(-aH --info=progress2 --partial)
if [[ "$DRY_RUN" -eq 1 ]]; then
    RSYNC_FLAGS+=(--dry-run)
fi

echo "Job ID:     ${SLURM_JOB_ID:-local}"
echo "Host:       $(hostname)"
echo "Source:     $SRC_ROOT"
echo "Dest:       $DEST_ROOT"
echo "Dry run:    $DRY_RUN"
echo "Extract:    $EXTRACT_TAR"
echo "Started:    $(date)"
echo
df -h "$DEST_ROOT" "$SRC_ROOT" 2>/dev/null || true
echo

[[ -d "$SRC_ROOT" ]] || { echo "Error: backup root missing: $SRC_ROOT" >&2; exit 1; }

fail=0
for mapping in "${RESTORE_DIRS[@]}"; do
    src_name="${mapping%%:*}"
    dest_name="${mapping##*:}"
    src="$SRC_ROOT/$src_name"
    dest="$DEST_ROOT/$dest_name"

    if [[ ! -d "$src" ]]; then
        echo "WARNING: skipping missing source: $src"
        continue
    fi

    echo "================================================================"
    echo "rsync  $src/  ->  $dest/"
    echo "size:  $(du -sh "$src" 2>/dev/null | awk '{print $1}')"
    echo "================================================================"

    mkdir -p "$dest"
    if rsync "${RSYNC_FLAGS[@]}" "$src/" "$dest/"; then
        echo "OK: $dest_name"
    else
        echo "FAILED: rsync $src_name (exit $?)" >&2
        fail=1
    fi
    echo
done

TAR_SRC="$DEST_ROOT/human_reference_proteome_liganded_pdbs_260529/masif_search.tar.gz"
TAR_DEST="$DEST_ROOT/masif_search"

if [[ "$EXTRACT_TAR" -eq 1 && "$DRY_RUN" -eq 0 && -f "$TAR_SRC" ]]; then
    echo "================================================================"
    echo "Extracting $TAR_SRC -> $TAR_DEST/"
    echo "================================================================"
    mkdir -p "$TAR_DEST"
    # Archive already has a top-level masif_search/ prefix.
    first_entry="$(tar -tzf "$TAR_SRC" | head -n 1 || true)"
    if [[ "$first_entry" == masif_search/* || "$first_entry" == "masif_search/" ]]; then
        tar -xzf "$TAR_SRC" -C "$DEST_ROOT"
    else
        tar -xzf "$TAR_SRC" -C "$TAR_DEST"
    fi
    echo "OK: extracted to $TAR_DEST"
    echo
elif [[ "$EXTRACT_TAR" -eq 1 && "$DRY_RUN" -eq 1 ]]; then
    echo "[dry-run] would extract $TAR_SRC -> $TAR_DEST/"
    echo
fi

echo "================================================================"
echo "Destination summary"
echo "================================================================"
du -sh "$DEST_ROOT"/* 2>/dev/null || true
echo
echo "Finished:  $(date)"
echo "Exit code: $fail"

if [[ "$fail" -ne 0 ]]; then
    echo "One or more rsyncs failed; re-run the same sbatch (rsync --partial resumes)." >&2
fi

echo
echo "Note: data/input, EvoEF2, and the singularity image are NOT in this"
echo "work backup. Reinstall / recreate those separately if needed."

exit "$fail"
