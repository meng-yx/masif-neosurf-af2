#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 6G
#SBATCH --time 1:00:00
#SBATCH --array=0-120
#SBATCH --output=exelogs/make_transformation_12A-%A/slurm-%A_%a.out
#SBATCH --error=exelogs/make_transformation_12A-%A/slurm-%A_%a.out

echo STARTING AT `date`

# Not needed because we are running inside singularity container
# source /work/upcorreia/bin/load_masif_environment.sh 
./make_transformations_12A.sh ${SLURM_ARRAY_TASK_ID}
echo FINISHED at `date`
