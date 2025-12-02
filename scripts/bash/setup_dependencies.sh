#! /bin/bash

# Copy singularity image to the repo root
repo_root=$(git rev-parse --show-toplevel)
cd $repo_root
echo "Copying singularity image to the repo root"
cp /work/upthomae/Meng/MaSIF/IMAGE/masif-neosurf_v0.1.sif .

# Create MaSIF conda environment if it doesn't exist
if [ ! -d "MaSIF" ]; then
    echo "Creating MaSIF conda environment"
    conda create -n MaSIF -f environment.yml
    conda activate MaSIF
else 
    echo "MaSIF conda environment already exists"
fi