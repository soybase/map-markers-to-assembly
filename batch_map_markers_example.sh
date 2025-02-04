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
source activate ds-curate

PATH=$PWD/bin:$PATH

# Map a marker set between two assemblies:
map-markers.sh -c config/gnm1_to_gnm2_SoySSR.conf

# map a marker set from one assembly into several others
for gnm in 1 2 4 6; do 
  CONFIG=gnm1_to_gnm${gnm}_SoySNP50K.conf
  echo "Working on $CONFIG"
  map-markers.sh -c config/$CONFIG
  echo
done

date   # print timestamp

