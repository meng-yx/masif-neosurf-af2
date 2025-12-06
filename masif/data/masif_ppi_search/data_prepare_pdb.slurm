#!/bin/bash
#SBATCH --account=lpdi
#SBATCH --nodes=1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7000
#SBATCH --time=01:00:00
#SBATCH --array=1
#SBATCH --output=exelogs/_masif_precompute.%A_%a.out
#SBATCH --error=exelogs/_masif_precompute.%A_%a.out

# Dynamically get repo root
masif_neosurf_root=$(git rev-parse --show-toplevel)
SH_SCRIPT="$masif_neosurf_root/masif/data/masif_ppi_search/data_prepare_one_pdb.sh"


i=1
while read p; do
    FIELD1=$(echo $p | cut -d" " -f1)
    PDBID=$(echo $FIELD1 | cut -d"_" -f1)
    CHAIN1=$(echo $FIELD1 | cut -d"_" -f2)
    CHAIN2=$(echo $FIELD1 | cut -d"_" -f3)

    # Run the .sh script with absolute path
    $SH_SCRIPT $PDBID\_$CHAIN1\_$CHAIN2
done < /home/ramanoel/masif-neosurf-af2/masif/data/masif_ppi_search/lists/small_list.txt
