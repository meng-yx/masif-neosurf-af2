#!/bin/bash -l

repo_root=$(git rev-parse --show-toplevel)
cd $repo_root

#-----------------------------------------------------
# Installation of EvoEF2
#-----------------------------------------------------

# Check if EvoEF2 already exists
if [ -d "EvoEF2" ]; then
    echo "EvoEF2 already exists"
else
    # Clone from github repo
    echo "Cloning EvoEF2 from github repo"
    git clone https://github.com/tommyhuangthu/EvoEF2.git

    # build EvoEF2.cpp to executable
    echo "Building EvoEF2.cpp to executable"
    cd EvoEF2
    chmod +x ./build.sh
    ./build.sh

    # Test if EvoEF2 executable is built successfully
    ./EvoEF2 --help

    # Got back to the root directory
    cd ..

    echo "EvoEF2 installation completed"
fi


