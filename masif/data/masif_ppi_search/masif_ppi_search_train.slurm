#!/bin/bash -l
#SBATCH --account=lpdi
#SBATCH --nodes=1
#SBATCH --time=40:0:0

#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard
#SBATCH --mem 7000
#SBATCH --output=exelogs/gpu_monet_seeder.%A_%a.out
#SBATCH --error=exelogs/gpu_monet_seeder.%A_%a.err

### #SBATCH --qos=gpu
module purge
module load gcc cuda cudnn mvapich2 openblas
deactivate
source ~/lpdi_fs/masif/tensorflow-1.12/bin/activate
./train_nn.sh nn_models.sc05.all_feat.custom_params
