#!/bin/bash
#SBATCH --time=03:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20  
#SBATCH --partition=short    # standard node(s)
#SBATCH --job-name="panfam"
#SBATCH --mail-user=steven.cannon@usda.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

set -o errexit
set -o nounset
set -o xtrace

date   # print timestamp

module load samtools bedtools blast+

export PATH=bin:$PATH
printf "Testing path: marker_gff_to_bed.pl is "; 
  which marker_gff_to_bed.pl
echo

marker_locs=markers/Wm82.gnm1.mrk.SoySNP50K/glyma.Wm82.gnm1.mrk.SoySNP50K.gff3.gz
genome_from=genomes/Wm82.gnm1.FCtY/glyma.Wm82.gnm1.FCtY.genome_main.fna
genome_to=genomes/Wm82.gnm2.DTC4/glyma.Wm82.gnm2.DTC4.genome_main.fna
marker_to_file=markers/Wm82.gnm2.mrk.SoySNP50K/glyma.Wm82.gnm2.mrk.SoySNP50K.tsv
threads=20

./map-markers-to-assembly.sh $marker_locs $genome_from $genome_to $marker_to_file $threads

date   # print timestamp


