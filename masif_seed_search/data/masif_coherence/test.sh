#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=00:01:00
#SBATCH --partition=h100
#SBATCH --gpus-per-node=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/test.%A.out
#SBATCH --error=logs/test.%A.out

conda activate masif-coherence
python -c "import torch; print(torch.__version__); print('cuda:', torch.cuda.is_available()); print(torch.version.cuda)"