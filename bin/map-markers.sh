#!/usr/bin/env bash

version="2025-01-08"

set -o errexit -o errtrace -o nounset -o pipefail -o posix

trap 'echo ${0##*/}:${LINENO} ERROR executing command: ${BASH_COMMAND}' ERR

# to help assign a heredoc value to a variable. The internal line return is intentional.
define(){ o=; while IFS=$'\n' read -r a; do o="$o$a"'
'; done; eval "$1=\$o"; }

########## Main program

scriptname='map-markers.sh'

define HELP_DOC <<'EOS'
NAME
  map-markers-to-assembly.sh  -- Given a file of marker locations in one genome, report the
    locations of those markers in a second genome.

SYNOPSIS
  map-markers.sh  -c CONFIG_FILE

  Required:
           -c (path to the config file)

  Options: -h help

  Specify paths to the "from" and "to" genome assemblies and the (gff3) marker file.
  The two genome files will be uncompressed in the work directory, as part of this script. The primary intended
  use of the script is on the same file system as the data files, in which case the starting files can be copied 
  into the work directory. If they are coming from a remote remote location, one solution would be to pull the 
  files locally into a data directory using scp or equivalent, and then give the paths to those files.

  VARIABLES set in config file:
    marker_from - Full filepath to file with marker names and locations on first Genome; in gff3 format, compressed

    genome_from - Full filepath to first genome assembly, corresponding with the coordinates in the marker_from file; compressed
    genome_to   - Full filepath to the second genome assembly, to which the markers will be projected

    marker_to   - Name for new marker file (gff3); will be written to work_dir/marker_to/

    identity    - Minimum percent identity in range 0..100 for blastn qcovhsp [90]
    sample_len  -  Maximum length of sequence variant to report, as a sample, in the GFF 9th column [10]
    max_len     - Maximum variant length for which to report a GFF line [200]

    work_dir    - work directory; default work_dir

AUTHOR
    Steven Cannon <steven.cannon@usda.gov>

EOS

if [ "$#" -eq 0 ]; then
  echo >&2 "$HELP_DOC" && exit 0;
fi

NPROC=$( ( command -v nproc > /dev/null && nproc ) || getconf _NPROCESSORS_ONLN)
CONFIG="null"

export NPROC=${NPROC:-1}

##########
# Command-line interpreter

while getopts "c:h" opt
do
  case $opt in
    c) CONFIG=$OPTARG; echo "Config: $CONFIG" ;;
    h) echo >&2 "$HELP_DOC" && exit 0 ;;
    *) echo >&2 echo "$HELP_DOC" && exit 1 ;;
  esac
done

if [ "$CONFIG" == "null" ]; then
  printf "\nPlease provide the path to a config file: -c CONFIG\n" >&2
  printf "\nRun \"%s -h\" for help.\n\n" "$scriptname" >&2
  exit 1;
else
  export CONF=${CONFIG}
fi

# Add shell variables from config file
# shellcheck source=/dev/null
. "${CONF}"

# Check for existence of third-party executables
missing_req=0
dependencies='blastn samtools bedtools'
for program in $dependencies; do
  if ! type "$program" &> /dev/null; then
    echo "Warning: executable $program is not on your PATH."
    missing_req=$((missing_req+1))
  fi
done
if [ "$missing_req" -gt 0 ]; then 
  printf "\nPlease add the programs above to your environment and try again.\n\n"
  exit 1; 
fi

# Check that the bin directory is in the PATH
if ! type marker_gff_to_bed_and_var.pl &> /dev/null; then
  printf "\nPlease add the bin directory to your PATH and try again. Try the following:\n"
  printf "\n  PATH=%s/bin:\%s\n\n" "$PWD" "$PATH"
  exit 1; 
fi

##########

cat_or_zcat() {
  case ${1} in
    *.gz) gzip -dc "$@" ;;
       *) cat "$@" ;;
  esac
}

##########

echo
echo "== Copy files into work directory. Uncompress and index where needed."

WD=`realpath $work_dir`
mkdir -p "${WD}"
cd "${WD}" || exit

mkdir -p marker_from
mkdir -p genome_from
mkdir -p genome_to
mkdir -p marker_to

# Copy files into working directory. Uncompress and index the "from" genome.
MRK_BASE=`basename $marker_from .gz`
GNM_FROM_BASE=`basename $genome_from .gz`
GNM_TO_BASE=`basename $genome_to .gz`

if [ ! -f  marker_from/$MRK_BASE.gz ]; then
  cp $marker_from marker_from/$MRK_BASE.gz || exit
fi

if [ ! -f genome_from/$GNM_FROM_BASE.gz ] || [ ! -f genome_from/$GNM_FROM_BASE ]; then
  cp $genome_from genome_from/$GNM_FROM_BASE.gz || exit

  if [ ! -f genome_from/$GNM_FROM_BASE ]; then
    gunzip genome_from/$GNM_FROM_BASE.gz
  fi
  
  if [ ! -f genome_from/$GNM_FROM_BASE.fai ]; then
    samtools faidx genome_from/$GNM_FROM_BASE
  fi
fi

if [ ! -f genome_to/$GNM_TO_BASE.gz ] || [ ! -f genome_to/$GNM_TO_BASE ]; then
  cp $genome_to genome_to/$GNM_TO_BASE.gz || exit

  if [ ! -f genome_to/$GNM_TO_BASE ]; then
    gunzip genome_to/$GNM_TO_BASE.gz
  fi
fi

MRK_FR_BARE="${MRK_BASE%.*}"
if [ ! -f marker_to/$MRK_FR_BARE.1kflank.fna ]; then
  echo
  echo "== Put the marker information into four-column BED format, with 1000 bases on each side of the SNP."
  echo "== Need to adjust from GFF 1-based, closed [start, end] to BED 0-based, half-open [start-1, end)."
  echo "== Special-casing near the molecule start with script marker_gff_to_bed_and_var.pl."
    marker_gff_to_bed_and_var.pl -marker $marker_from \
                                 -genome genome_from/$GNM_FROM_BASE \
                                 -out "$WD/marker_from/$MRK_FR_BARE"
  
  echo
  echo "== Extract the sequences from the FROM genome"

  bedtools getfasta -fi genome_from/$GNM_FROM_BASE \
                    -bed marker_from/$MRK_FR_BARE.UD.bed \
                    -name \
                    -fo marker_to/$MRK_FR_BARE.1kflank.fna
  
  echo
  echo "== Strip positional information, added by getfasta, from the retrieved sequence"
  perl -pi -e 's/>(\S+)::.+/>$1/' marker_to/$MRK_FR_BARE.1kflank.fna
else 
  echo "Skipping creation of BED file and extraction of flanking sequence, since it exists."
fi

echo
echo "== Make BLAST output directories and index files"
mkdir -p blastdb blastout

echo "genome_from: genome_to/$GNM_TO_BASE"
if [ ! -f blastdb/$GNM_TO_BASE.nin ]; then
  cat_or_zcat genome_to/$GNM_TO_BASE | 
    makeblastdb -in - -dbtype nucl \
                -title $GNM_TO_BASE -hash_index -out blastdb/$GNM_TO_BASE
fi

if [ ! -f blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bln ]; then
  echo
  echo "== Run BLAST"
  blastn -db blastdb/$GNM_TO_BASE \
         -query marker_to/$MRK_FR_BARE.1kflank.fna \
         -num_threads $NPROC -evalue 1e-10 -perc_identity 99 \
         -outfmt "6 std qlen qcovs" \
         -out blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bln
else 
  echo "== Skipping BLAST, as the output file exists."
fi

echo
echo "== Filter BLAST output and write new marker file (as a tsv file)"
cat blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bln | top_line.awk | 
  marker_blast_to_gff.pl -genome genome_to/$GNM_TO_BASE \
                         -gff_source $gff_source \
                         -gff_type $gff_type \
                         -gff_ID_prefix $gff_ID_prefix \
			 -max_len $max_len \
                         -out "$WD/marker_to/$marker_to"

echo
echo "== Sort GFF"
sort_gff.pl "$WD/marker_to/$marker_to.gff3" > "$WD/marker_to/tmp.gff"
mv "$WD/marker_to/tmp.gff" "$WD/marker_to/$marker_to.gff3"


echo
echo "== Compare the initial and mapped markers and report"
echo "==   Marker list 1: $work_dir/marker_to/lis.$MRK_FR_BARE"
echo "==   Marker list 2: $work_dir/marker_to/lis.$marker_to"
echo "==   Marker report: $work_dir/marker_to/report.${MRK_FR_BARE}--${marker_to}.tsv"

# Extract ID and allele from the "from" bed file, and print "+" orientation for all
cat "$WD/marker_from/$MRK_FR_BARE.bed" | awk -v OFS="\t" '{print $4, "+", $5}' | sort > "$WD/marker_to/lis.$MRK_FR_BARE"

# Next: Extract ID, allele, and orientation from the "to" bed file
cut -f4,6,7 "$WD/marker_to/$marker_to.bed" | sort > "$WD/marker_to/lis.$marker_to"

# Fields in joined result: markerID, orient, allele_gnm1, orient, allele_gnm2
#                              0        1  2  3  4
#                          ss715578401  +  G  +  G
#                          ss715578490  +  G  -  C
join -a1 "$WD/marker_to/lis.$MRK_FR_BARE" "$WD/marker_to/lis.$marker_to" |
  perl -F"\s" -lane 'BEGIN{ print join("\t", "#markerID", "compare", "len1", "len2", "orient", "var1", "var2") };
                     if (scalar(@F)==3){ # marker not in target genome
                       print join( "\t", $F[0], "NULL", length($F[2]), 0, ".", $F[2], "NULL");
                     }
                     else {
                       if ($F[2] eq $F[4]){ # allele is the same in both genomes
                         print join( "\t", $F[0], "same", length($F[2]), length($F[4]), $F[3], $F[2], $F[4]);
                       }
                       else {
                         print join( "\t", $F[0], "NOT",  length($F[2]), length($F[4]), $F[3], $F[2], $F[4]);
                       }
                     }
                    '  > "$WD/marker_to/report.${MRK_FR_BARE}--${marker_to}.tsv"

# TO DO: get allele from bed file rather than gff, or add orientation to bed file

echo
echo "== Generate report of marker orientations"
cat "$WD/marker_to/$marker_to.gff3" | sort -k1,1 -k4n,4n |
   awk -v ORS=" " '$1 == prev {print $8; prev=$1} 
                   NR!=1 && $1 != prev {print "\n\n" $1 "\n" $8 ; prev=$1} 
                   NR==1 {print $1 "\n" $8 ; prev=$1}
                   END{print "\n"}' > "$WD/marker_to/orient.${MRK_FR_BARE}--${marker_to}.txt"

echo
echo "== Mapped markers:"
echo "==   $work_dir/marker_to/$marker_to.bed"
echo "==   $work_dir/marker_to/$marker_to.gff3"
echo
echo "== Run completed. Look for results at $work_dir/marker_to/"
echo


exit 0
