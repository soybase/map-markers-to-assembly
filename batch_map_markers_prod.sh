#!/bin/bash
#SBATCH --time=23:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # 20 processor core(s) per node X 2 threads per core
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="map-markers"
#SBATCH --mail-user=steven.cannon@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

set -o errexit
set -o nounset
#set -o xtrace

date   # print timestamp

# If using conda environment for dependencies:
module load miniconda
source activate ds-curate

PATH=$PWD/bin:$PATH

map-markers.sh -c config/gnm1_to_Fisk_SoySNP50K.conf

# # Map SoySSR set
# for gnm in 1 2 4 6; do 
#   CONFIG=gnm1_to_gnm${gnm}_SoySSR.conf
#   echo "Working on $CONFIG"
#   map-markers.sh -c config/$CONFIG
#   echo
# done
# 
# # map SoySNP50K set
# for gnm in 1 2 4 6; do 
#   CONFIG=gnm1_to_gnm${gnm}_SoySNP50K.conf
#   echo "Working on $CONFIG"
#   map-markers.sh -c config/$CONFIG
#   echo
# done

date   # print timestamp

