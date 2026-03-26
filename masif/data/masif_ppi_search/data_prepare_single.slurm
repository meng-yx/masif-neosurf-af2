#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=bigmem
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 2
#SBATCH --mem 28000
#SBATCH --time 01:00:00
#SBATCH --output=exelogs/masif_precompute_single_%A.out
#SBATCH --error=exelogs/masif_precompute_single_%A.out

COMPLEX_ID="$1"
PINDER_PDB_DIR="$2"
if [ -z "$COMPLEX_ID" ] || [ -z "$PINDER_PDB_DIR" ]; then
    echo "Usage: sbatch data_prepare_single.slurm COMPLEX_ID /path/to/pinder_pdb_dir"
    exit 1
fi
./data_prepare_one.sh $COMPLEX_ID "$PINDER_PDB_DIR"
