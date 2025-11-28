#!/bin/bash
#SBATCH --nodes 1
#SBATCH --partition=standard
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 7000
#SBATCH --time 01:00:00
#SBATCH --array=1
#SBATCH --output=exelogs/compute_descriptors.%A_%a.out
#SBATCH --error=exelogs/compute_descriptors.%A_%a.out

# TODO: See if it runs faster on GPU

i=1
while read p; do
    if [ $(( i % 1001 )) == ${SLURM_ARRAY_TASK_ID} ]; then
        FIELD1=$(echo $p| cut -d" " -f1)
        PDBID=$(echo $FIELD1| cut -d"_" -f1)
        CHAIN1=$(echo $FIELD1| cut -d"_" -f2)
        CHAIN2=$(echo $FIELD1| cut -d"_" -f3)
        echo "Computing descriptors for $PDBID\_$CHAIN1\_$CHAIN2"
        ./compute_descriptors.sh $PDBID\_$CHAIN1\_$CHAIN2
    fi
    i=$((i+1))
done < lists/full_list.txt
