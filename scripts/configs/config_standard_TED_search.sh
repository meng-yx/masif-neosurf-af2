################################################################################
# Centralized configuration for pipeline scripts
# Source from repo root:
#   source scripts/config.sh
################################################################################

# ----------------- General settings -----------------

# MaSIF-search defaults
CUTOFF=6.0
IFACE_CUTOFF=0.1
DESC_DIST_CUTOFF=3
NN_SCORE_CUTOFF=0.9
RANSAC_ITER=100000
N_RETRY_ALIGNMENT=1
TARGET_SAMPLING_RATIO=0.2

# Neosurf-Neosurf settings (automatically set anchor residue to the largest HETAOM residue)
TARGET_AUTO_NEOSURF=1
SEED_AUTO_NEOSURF=0
RESUME=1

# Utility
PIP_PACKAGES="pandas tqdm"


################################################################################
# Paths (relative to masif-neosurf repo root)
################################################################################

IMAGE=masif-neosurf_v0.1.sif
BIND_MOUNTS="/work:/work,/scratch:/scratch"


# Seed database root
#DATABASE_ROOT=data/preprocess
DATABASE_ROOT=/work/upthomae/Meng/TED_human_domainome_MaSIF/output

# Target database root
TARGET_PREPROCESS_ROOT=data/preprocess

REDUCE_HET_DICT=/work/lpdi/users/schneuing/misc/reduce_wwPDB_het_dict.txt
