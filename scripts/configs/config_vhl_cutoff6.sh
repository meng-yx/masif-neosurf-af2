################################################################################
# Configuration preset: VHL-specific neosurf-neosurf search at CUTOFF = 6.
#
# Identical to scripts/config.sh EXCEPT CUTOFF is 6.0 (instead of 4.5), so both
# the target and seed neosurface patches are sampled further from the ligand
# (--target_site_cutoff / --seed_site_cutoff). Preprocessing is cutoff-
# independent, so this preset reuses the same DATABASE_ROOT / descriptors.
#
# Select it by exporting CONFIG_FILE before submitting, e.g.:
#   CONFIG_FILE=scripts/configs/config_vhl_cutoff6.sh sbatch scripts/slurm/search_array.sh ...
################################################################################

# ----------------- General settings -----------------

# MaSIF-search defaults
CUTOFF=6                       # <-- widened from 4.5
IFACE_CUTOFF=0.1
DESC_DIST_CUTOFF=3
NN_SCORE_CUTOFF=0.95
RANSAC_ITER=50000
N_RETRY_ALIGNMENT=1
TARGET_SAMPLING_RATIO=0.15

# Neosurf-Neosurf settings (automatically set anchor residue to the largest HETATM residue)
TARGET_AUTO_NEOSURF=1
SEED_AUTO_NEOSURF=1
RESUME=1

# Utility
PIP_PACKAGES="pandas tqdm"


################################################################################
# Paths (relative to masif-neosurf repo root)
################################################################################

IMAGE=masif-neosurf_v0.1.sif
BIND_MOUNTS="/work:/work,/scratch:/scratch"

# Seed database root (reused: surfaces/descriptors are cutoff-independent)
DATABASE_ROOT=data/preprocess

# Target database root
TARGET_PREPROCESS_ROOT=data/preprocess

REDUCE_HET_DICT=/work/lpdi/users/schneuing/misc/reduce_wwPDB_het_dict.txt
