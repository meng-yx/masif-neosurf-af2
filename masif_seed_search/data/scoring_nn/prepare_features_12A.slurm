#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 6G
#SBATCH --time 1:00:00
#SBATCH --array=1-120
#SBATCH --output=exelogs/prepare_features_12A-%A/slurm-%A_%a.out
#SBATCH --error=exelogs/prepare_features_12A-%A/slurm-%A_%a.out

echo STARTING AT `date`
echo ${SLURM_ARRAY_TASK_ID}

i=1
while read p; do
    if [ $(( i % 1001 )) == ${SLURM_ARRAY_TASK_ID} ]; then
        echo $p
        FIELD1=$(echo $p| cut -d" " -f1)
        PDBID=$(echo $FIELD1| cut -d"_" -f1)
        CHAIN1=$(echo $FIELD1| cut -d"_" -f2)
        CHAIN2=$(echo $FIELD1| cut -d"_" -f3)
        ./precompute_evaluation_features_12A.sh $PDBID\_$CHAIN1\_$CHAIN2
    fi
    i=$((i+1))
done < lists/full_list.txt

echo FINISHED at `date`
