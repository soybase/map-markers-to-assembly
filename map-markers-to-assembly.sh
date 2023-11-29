#!/usr/bin/env bash
# NAME
#   map-markers-to-assembly.sh  -- Given a file of marker locations in one genome, report the 
#     locations of those markers in a second genome.
#
# SYNOPSIS
#   map-markers-to-assembly.sh $marker_locs $genome_from $genome_to marker_to_file
#
# OPERANDS
#   Paths to four files (including the new marker-locations file to be created),
#   and number of threads to use in blast search
#     marker_locs    - Path to file of marker names and locations on Genome assembly 1; in gff3 format
#     genome_from    - Path to genome assembly 1, corresponding with the coordinates in the marker_locs file
#     genome_to      - Path to genome assembly 2, to which the markers will be projected
#     marker_to_file - Path for the marker-locations file to be created; may include a directory path
#     threads        - Number of threads to use in blast search
#
# AUTHOR
#     Steven Cannon <steven.cannon@usda.gov>

set -o errexit
set -o nounset

marker_locs=$1
genome_from=$2
genome_to=$3
marker_to_file=$4
threads=$5

##### Subroutines #####

cat_or_zcat() {
  case ${1} in
    *.gz) gzip -dc "$@" ;;
       *) cat "$@" ;;
  esac
}

##### Main #####

if [[ "$genome_from" =~ ".gz" ]]; then
  gunzip $genome_from
  FROM_GNM_BASE=`basename $genome_from .gz`
else
  FROM_GNM_BASE=`basename $genome_from`
fi
ext1="${FROM_GNM_BASE##*.}"
FROM_GNM_BASE=`basename $FROM_GNM_BASE .$ext1`
echo "== FROM_GNM_BASE: $FROM_GNM_BASE  .$ext1"

if [[ "$genome_to" =~ ".gz" ]]; then
  gunzip $genome_to
  TO_GNM_BASE=`basename $genome_to .gz`
else
  TO_GNM_BASE=`basename $genome_to`
fi
ext2="${TO_GNM_BASE##*.}"
TO_GNM_BASE=`basename $TO_GNM_BASE .$ext2`
echo "== TO_GNM_BASE: $TO_GNM_BASE  .$ext2"

if [[ "$marker_locs" =~ ".gz" ]]; then
  MRK_BASE=`basename $marker_locs .gz`
else
  MRK_BASE=`basename $marker_locs`
fi
ext3="${MRK_BASE##*.}"
MRK_BASE=`basename $MRK_BASE .$ext3`
echo "== MRK_BASE: $MRK_BASE  .$ext3"
echo

MRK_DIR=`dirname $marker_locs`
FROM_GNM_DIR=`dirname $genome_from`
TO_GNM_DIR=`dirname $genome_to`

echo "== Create an index on the uncompressed FROM genome"
echo
samtools faidx $FROM_GNM_DIR/$FROM_GNM_BASE.$ext1

echo "== Put the marker information into four-column BED format, with 1000 bases on each side of the SNP."
echo "== Need to adjust from GFF 1-based, closed [start, end] to BED 0-based, half-open [start-1, end)."
echo "== Special-casing near the molecule start with script marker_gff_to_bed.pl."
echo
cat_or_zcat $marker_locs |
  marker_gff_to_bed.pl -fai $FROM_GNM_DIR/$FROM_GNM_BASE.$ext1.fai \
    -out $MRK_DIR/$MRK_BASE.bed

echo "== Extract the sequences from the FROM genome"
echo
bedtools getfasta -fi $FROM_GNM_DIR/$FROM_GNM_BASE.$ext1 \
                  -bed $MRK_DIR/$MRK_BASE.bed \
                  -name \
                  -fo $MRK_DIR/$MRK_BASE.1kflank.fna

echo "== Strip positional information, added by getfasta, from the retrieved sequence"
perl -pi -e 's/>(\S+)::.+/>$1/' $MRK_DIR/$MRK_BASE.1kflank.fna

echo "== Make BLAST output directories and index files"
mkdir -p blastdb blastout  
makeblastdb -in $TO_GNM_DIR/$TO_GNM_BASE.$ext1 -dbtype nucl -hash_index -out blastdb/$TO_GNM_BASE

echo "== Run BLAST"
blastn -db blastdb/$TO_GNM_BASE \
       -query $MRK_DIR/$MRK_BASE.1kflank.fna \
       -num_threads $threads -evalue 1e-10 -perc_identity 99 \
       -outfmt "6 std qlen qcovs" \
       -out blastout/$MRK_BASE.x.$TO_GNM_BASE.bln 

MRK_TO_BASE=`basename $marker_to_file`
MRK_TO_DIR=`dirname $marker_to_file`

mkdir -p $MRK_TO_DIR

echo "== Filter BLAST output and write new marker file (as a tsv file)"
cat blastout/$MRK_BASE.x.$TO_GNM_BASE.bln | 
  awk '$4>=90 && $14>=90' | top_line.awk | sort -k2,2 -k9n,9n | 
  filter_marker_blast_data.awk > $MRK_TO_DIR/$MRK_TO_BASE

# TO DO: Reformat into gff3, with user-specified molecule prefix (col1) and source (col2).
