#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=96:00:00
#SBATCH --partition=h100
#SBATCH --qos=long
#SBATCH --gpus-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=exelogs/masif_ppi_search_train.%A.out
#SBATCH --error=exelogs/masif_ppi_search_train.%A.out

module purge
module load gcc cuda cudnn openblas
#conda deactivate

./train_nn.sh nn_models.sc05.all_feat.custom_params
