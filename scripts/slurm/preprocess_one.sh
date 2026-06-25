#!/bin/bash
#SBATCH --job-name=preprocess
#SBATCH --output=logs/preprocess_one-%A.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=28000


# PDB_PATH=$1
# TARGET=$2
# LIGAND=$3
# LIGAND_PATH=$4
# OUTPUT_DIR=$5

PDB_PATH=/scratch/ymeng/Neosurf_Neosurf/data/input/8VLB-3JF_A.pdb
TARGET=8VLB-3JF_A
LIGAND=3JF_A
LIGAND_PATH=/scratch/ymeng/Neosurf_Neosurf/data/input/8VLB_A_3JF.sdf
OUTPUT_DIR=/scratch/ymeng/Neosurf_Neosurf/data/preprocess

bash scripts/bash/preprocess.sh $PDB_PATH $TARGET $LIGAND $LIGAND_PATH $OUTPUT_DIR

