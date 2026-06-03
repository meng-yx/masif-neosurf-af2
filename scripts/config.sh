################################################################################
# Centralized configuration for pipeline scripts
# Source from repo root:
#   source scripts/config.sh
################################################################################

# ----------------- General settings -----------------

# MaSIF-search defaults
CUTOFF=6.0
IFACE_CUTOFF=0.3
DESC_DIST_CUTOFF=2.0
NN_SCORE_CUTOFF=0.9
RANSAC_ITER=100000
N_RETRY_ALIGNMENT=1
TARGET_SAMPLING_RATIO=0.2

# Utility
PIP_PACKAGES="pandas tqdm"

# Array job settings
ARRAY_JOBS_N=500
ARRAY_JOBS_RANGE="1-${ARRAY_JOBS_N}"

################################################################################
# Paths (relative to masif-neosurf repo root)
################################################################################

IMAGE=masif-neosurf_v0.1.sif
BIND_MOUNTS="/work:/work,/scratch:/scratch"

DATABASE_ROOT=data/preprocess

REDUCE_HET_DICT=/work/lpdi/users/schneuing/misc/reduce_wwPDB_het_dict.txt
