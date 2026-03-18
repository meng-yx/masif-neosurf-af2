#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=40:0:0
#SBATCH --cpus-per-task=4
#SBATCH --mem 24G
#SBATCH --output=exelogs/masif_ppi_search_train.%A_%a.out
#SBATCH --error=exelogs/masif_ppi_search_train.%A_%a.out

module purge
module load gcc cuda cudnn mvapich2 openblas
deactivate
source ~/lpdi_fs/masif/tensorflow-1.12/bin/activate
./train_nn.sh nn_models.sc05.all_feat.custom_params
