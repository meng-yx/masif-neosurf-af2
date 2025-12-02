#!/bin/bash

# Step 1: Generate Alignment Training Data (Parallelized)

# For 12 Angstrom patches
sbatch make_transformations_12A.slurm

# Step 2: Precompute Evaluation Features (Parallelized)

# For 12 Angstrom patches
sbatch prepare_features_12A.slurm

# Step 3: Train Neural Network
sbatch train_nn.slurm