#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 7000
#SBATCH --time 01:00:00
#SBATCH --output=exelogs/masif_precompute_%A/_masif_precompute.%A_%a.out
#SBATCH --error=exelogs/masif_precompute_%A/_masif_precompute.%A_%a.out

echo "Starting at $(date)"
PINDER_PDB_DIR="$1"
if [ -z "$PINDER_PDB_DIR" ]; then
    echo "Usage: sbatch data_prepare.slurm /path/to/pinder_pdb_dir"
    exit 1
fi

# # Check if the previous job has been OOM killed by searching for "OOM"
# last_job_output="/scratch/ymeng/masif-neosurf-pinder/masif/data/masif_ppi_search/exelogs/masif_precompute_52632200/_masif_precompute.52632200_${SLURM_ARRAY_TASK_ID}.out"
# if [ -f "$last_job_output" ]; then
#     if ! grep -q "OOM" "$last_job_output"; then
#         echo "Previous job did NOT run out of memory (OOM not found in $last_job_output), terminating this job."
#         exit 1
#     fi
# fi

# echo "Previous job was OOM killed, continuing with next job."

i=1
while read p; do
    if [ $(( i % 1001 )) == ${SLURM_ARRAY_TASK_ID} ]; then
        FIELD1=$(echo $p| cut -d" " -f1)
        PDBID=$(echo $FIELD1| cut -d"_" -f1)
        CHAIN1=$(echo $FIELD1| cut -d"_" -f2)
        CHAIN2=$(echo $FIELD1| cut -d"_" -f3)
        ./data_prepare_one.sh $PDBID\_$CHAIN1\_$CHAIN2 "$PINDER_PDB_DIR"
    fi
    i=$((i+1))
done < lists/full_list.txt
echo "Finished at $(date)"