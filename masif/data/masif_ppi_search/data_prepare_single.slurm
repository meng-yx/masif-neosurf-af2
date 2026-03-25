#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 7000
#SBATCH --time 01:00:00
#SBATCH --array=1
#SBATCH --output=exelogs/_masif_precompute.%A_%a.out
#SBATCH --error=exelogs/_masif_precompute.%A_%a.out

PDBID="1F3X"
CHAIN1="G"
CHAIN2="H"
PINDER_PDB_DIR="$1"
if [ -z "$PINDER_PDB_DIR" ]; then
    echo "Usage: sbatch data_prepare_single.slurm /path/to/pinder_pdb_dir"
    exit 1
fi
./data_prepare_one.sh $PDBID\_$CHAIN1\_$CHAIN2 "$PINDER_PDB_DIR"
