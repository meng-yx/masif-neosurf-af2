#!/bin/bash
jid1=$(sbatch data_prepare.slurm | awk '{print $4}')
jid2=$(sbatch --dependency=afterok:$jid1 cache_nn.slurm | awk '{print $4}')
jid3=$(sbatch --dependency=afterok:$jid2 masif_ppi_search_train.slurm | awk '{print $4}')
jid4=$(sbatch --dependency=afterok:$jid3 compute_descriptors.slurm | awk '{print $4}')
