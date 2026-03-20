#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --partition=h100
#SBATCH --qos=normal
#SBATCH --gpus-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=exelogs/masif_ppi_search_train.%A_%a.out
#SBATCH --error=exelogs/masif_ppi_search_train.%A_%a.out

module purge
module load gcc cuda cudnn openblas
#conda deactivate

./train_nn.sh nn_models.sc05.all_feat.custom_params
