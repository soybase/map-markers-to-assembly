#!/bin/bash
#SBATCH   # Your ...
#SBATCH   # sbatch ...
#SBATCH   # commands ...
#SBATCH   # here ...
#SBATCH   # ...

set -o errexit
set -o nounset
#set -o xtrace

date   # print timestamp

# If using conda environment for dependencies:
module load conda
source activate map-markers

PATH=$PWD/bin:$PATH

# Map a marker set between two assemblies:
map-markers.sh -c config/gnm1_to_gnm2_SoySSR.conf

date   # print timestamp

