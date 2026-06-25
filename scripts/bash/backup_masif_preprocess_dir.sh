#!/usr/bin/env bash
#
# Back up only the files required to run masif_search.py.
#
# Keeps the union of seed- and target-required artifacts:
#   - PDBs, benchmark surfaces, MaSIF-site predictions, descriptors
#   - 04b precomputation patch indices only (not 04a or other 04b arrays)
#   - pred_surfaces PLYs (target-only at search time; kept for all proteins
#     because the target set can change)
#
# Modes:
#   rsync    Copy files into DEST (use before deleting the original tree).
#   symlink  Mirror layout with symlinks into DEST (quick test while SOURCE remains).
#
# Usage (from repo root or any cwd):
#   scripts/bash/backup_masif_preprocess_dir.sh rsync  data/preprocess  data/preprocess_search
#   scripts/bash/backup_masif_preprocess_dir.sh symlink data/preprocess  data/preprocess_search
#   scripts/bash/backup_masif_preprocess_dir.sh rsync --dry-run data/preprocess /tmp/test
#
set -euo pipefail

MODE=""
SOURCE=""
DEST=""
DRY_RUN=0

usage() {
    cat <<'EOF'
Usage: backup_masif_preprocess_dir.sh MODE SOURCE DEST [--dry-run]

MODE:
  rsync     Copy search-required files into DEST.
  symlink   Create DEST tree with symlinks to SOURCE (for testing).

SOURCE:    MaSIF preprocess root (e.g. data/preprocess).
DEST:      Output directory (created if missing).

Only files needed by masif_search.py are included. Everything under
04a-precomputation_9A/ and non-list_indices arrays under 04b are skipped.
EOF
}

die() {
    echo "Error: $*" >&2
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        rsync|symlink) MODE="$1" ;;
        --dry-run) DRY_RUN=1 ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            if [[ -z "$SOURCE" ]]; then
                SOURCE="$1"
            elif [[ -z "$DEST" ]]; then
                DEST="$1"
            else
                die "unexpected argument: $1"
            fi
            ;;
    esac
    shift
done

[[ -n "$MODE" ]] || { usage; exit 1; }
[[ -n "$SOURCE" ]] || die "SOURCE is required"
[[ -n "$DEST" ]] || die "DEST is required"

SOURCE="$(readlink -f "$SOURCE")"
DEST="$(readlink -f -m "$DEST")"

[[ "$SOURCE" != "$DEST" ]] || die "SOURCE and DEST must differ"
[[ -d "$SOURCE" ]] || die "SOURCE is not a directory: $SOURCE"

# Minimal sanity check: expected preprocess layout.
for required in \
    "data_preparation/01-benchmark_pdbs" \
    "data_preparation/01-benchmark_surfaces" \
    "data_preparation/04b-precomputation_12A/precomputation" \
    "descriptors/sc05/all_feat" \
    "output/all_feat_3l/pred_data" \
    "output/all_feat_3l/pred_surfaces"
do
    [[ -d "$SOURCE/$required" ]] || die "SOURCE does not look like a preprocess root (missing $required)"
done

collect_relative_paths() {
    local root="$1"
    (
        cd "$root"
        find data_preparation/01-benchmark_pdbs -type f -name '*.pdb' 2>/dev/null
        find data_preparation/01-benchmark_surfaces -type f -name '*.ply' 2>/dev/null
        find output/all_feat_3l/pred_data -type f -name 'pred_*.npy' 2>/dev/null
        find output/all_feat_3l/pred_surfaces -type f -name '*.ply' 2>/dev/null
        find descriptors/sc05/all_feat -type f \( -name 'p*_desc_straight.npy' -o -name 'p*_desc_flipped.npy' \) 2>/dev/null
        find data_preparation/04b-precomputation_12A/precomputation -type f -name '*_list_indices.npy' 2>/dev/null
    ) | LC_ALL=C sort -u
}

FILE_LIST="$(mktemp)"
trap 'rm -f "$FILE_LIST"' EXIT
collect_relative_paths "$SOURCE" >"$FILE_LIST"

N_FILES=$(wc -l <"$FILE_LIST" | tr -d ' ')
[[ "$N_FILES" -gt 0 ]] || die "no search-required files found under $SOURCE"

echo "Mode:       $MODE"
echo "Source:     $SOURCE"
echo "Dest:       $DEST"
echo "Files:      $N_FILES"
echo "Dry run:    $DRY_RUN"
echo

backup_rsync() {
    if [[ "$DRY_RUN" -eq 1 ]]; then
        rsync -av --dry-run --files-from="$FILE_LIST" "$SOURCE/" "$DEST/"
        return
    fi

    mkdir -p "$DEST"
    rsync -a --info=progress2 --files-from="$FILE_LIST" "$SOURCE/" "$DEST/"
}

backup_symlink() {
    local rel src_path dest_path dest_dir

    while IFS= read -r rel; do
        [[ -n "$rel" ]] || continue
        src_path="$SOURCE/$rel"
        dest_path="$DEST/$rel"
        dest_dir="$(dirname "$dest_path")"

        if [[ "$DRY_RUN" -eq 1 ]]; then
            echo "ln -sf $src_path $dest_path"
            continue
        fi

        mkdir -p "$dest_dir"
        ln -sf "$src_path" "$dest_path"
    done <"$FILE_LIST"
}

case "$MODE" in
    rsync) backup_rsync ;;
    symlink) backup_symlink ;;
esac

if [[ "$DRY_RUN" -eq 0 ]]; then
    echo
    echo "Done. $N_FILES paths written under $DEST"
    if [[ "$MODE" == "symlink" ]]; then
        echo "Symlink tree points at $SOURCE — keep SOURCE until you rsync-copy and verify search."
    else
        echo "Point masif_search --seed_dir / --target_dir at $DEST and verify before removing $SOURCE."
    fi
fi
