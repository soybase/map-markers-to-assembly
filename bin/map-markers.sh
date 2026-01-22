#!/usr/bin/env bash

version="2026-01-21"

 set -xtrace  # uncomment for debugging
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

    qcov_identity    - Minimum percent identity in range 0..100 for blastn qcovhsp [90]
    perc_identity    - Minimum percent identity in range 0..100 for blastn sequence match [99]
    sample_len  - Maximum length of sequence variant to report, as a sample, in the GFF 9th column [10]
    max_len     - Maximum variant length for which to report a GFF line [200]
    engine      - blast or burst [blast]
                    BLAST should work well in essentially every situation. The reason to consider BURST is that it is
                    much faster for very large marker sets. The downsides to BURST are slightly lower sensitivity, 
                    and the target-genome index files are about 20x larger than for BLAST; and their creation takes 
                    a large amount of memory (500 GB for a typical genome).

    work_dir    - work directory; default work_dir

AUTHOR
    Steven Cannon <steven.cannon@usda.gov>

EOS

if [ "$#" -eq 0 ]; then
  echo >&2 "$HELP_DOC" && exit 0;
fi

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
# Some helper functions ...
cat_or_zcat() {
  case ${1} in
    *.gz) gzip -dc "$@" ;;
       *) cat "$@" ;;
  esac
}

# print seqid length
seqlen() {
  fastafile=${1}
  awk '/^>/ {if (len) print len; len = 0; printf("%s\t", substr($1,2)); next}
       { len+=length }
       END {if (len) print len}' $fastafile
}
##########

echo
echo "Run of the map-markers-to-genome workflow, version $version"
echo "== Copy files into work directory. Uncompress and index where needed."

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

# Add shell variables from config file. Defaults, overridden by config file
marker_from=""; genome_from=""; genome_to=""; marker_to=""; gff_source=""; gff_ID_prefix=""; 
engine="burst"; gff_type="genetic_marker"; perc_identity="95"; gff_prefix_regex='^[^.]+\.[^.]+\.[^.]+\.'; 
evalue="1e-10"; qcov_identity="90"; sample_len="10"; max_len="25"; work_dir="work_dir"; min_flank="100";
# shellcheck source=/dev/null
. "${CONF}"

WD=$(realpath "$work_dir")
mkdir -p "${WD}"

mkdir -p "${WD}/marker_from"
mkdir -p "${WD}/genome_from"
mkdir -p "${WD}/genome_to"
mkdir -p "${WD}/marker_to"

# Copy files into working directory. Uncompress and index the "from" genome. Strip ".gz"
MRK_BASE=$(basename "$marker_from" .gz)
GNM_FROM_BASE=$(basename "$genome_from" .gz)
GNM_TO_BASE=$(basename "$genome_to" .gz)

# Strip suffix such as .fa, .fasta, .fna
GNM_FROM_BASE="${GNM_FROM_BASE%.*}"
GNM_TO_BASE="${GNM_TO_BASE%.*}"

if [ ! -f  "${WD}/marker_from/$MRK_BASE".gz ]; then
  cp "$marker_from" "${WD}/marker_from/$MRK_BASE".gz || exit
fi

if [ ! -f "${WD}/genome_from/$GNM_FROM_BASE".gz ] || [ ! -f "${WD}/genome_from/$GNM_FROM_BASE" ]; then
  cp "$genome_from" "${WD}/genome_from/$GNM_FROM_BASE".gz || exit

  if [ ! -f "${WD}/genome_from/$GNM_FROM_BASE" ]; then
    gunzip "${WD}/genome_from/$GNM_FROM_BASE".gz
  fi
  
  if [ ! -f "${WD}/genome_from/$GNM_FROM_BASE".fai ]; then
    samtools faidx "${WD}/genome_from/$GNM_FROM_BASE"
  fi
fi

if [ ! -f "${WD}/genome_to/$GNM_TO_BASE".gz ] || [ ! -f "${WD}/genome_to/$GNM_TO_BASE" ]; then
  cp "$genome_to" "${WD}/genome_to/$GNM_TO_BASE".gz || exit

  if [ ! -f "${WD}/genome_to/$GNM_TO_BASE" ]; then
    gunzip "${WD}/genome_to/$GNM_TO_BASE".gz
  fi
fi

MRK_FR_BARE="${MRK_BASE%.*}"
if [ ! -f "${WD}/marker_from/$MRK_FR_BARE".1kflank.fna ]; then
  echo
  echo "== Put the marker information into four-column BED format, with 1000 bases on each side of the SNP."
  echo "== Need to adjust from GFF 1-based, closed [start, end] to BED 0-based, half-open [start-1, end)."
  echo "== Special-casing near the molecule start with script marker_gff_to_bed_and_var.pl."
    marker_gff_to_bed_and_var.pl -marker "$marker_from" \
                                 -min_flank "$min_flank" \
                                 -genome "$WD/genome_from/$GNM_FROM_BASE" \
                                 -out "$WD/marker_from/$MRK_FR_BARE"
  
  echo
  echo "== Extract the sequences from the FROM genome"

  bedtools getfasta -fi "$WD/genome_from/$GNM_FROM_BASE" \
                    -bed "$WD/marker_from/$MRK_FR_BARE".UD.bed \
                    -name | 
                      filter_fasta_for_Ns.awk -v min_nonN="$min_flank" > "$WD/marker_from/$MRK_FR_BARE.1kflank.fna"
  
  echo
  echo "== Strip positional information, added by getfasta, from the retrieved sequence"
  perl -pi -e 's/>(\S+)::.+/>$1/' "$WD/marker_from/$MRK_FR_BARE".1kflank.fna
else 
  echo "Skipping creation of BED file and extraction of flanking sequence, since it exists."
fi

echo
echo "== Make sequence-search output directories and index files for search engine $engine"
mkdir -p "$WD/blastdb" "$WD/blastout"
  echo "genome_to: genome_to/$GNM_TO_BASE"

if [[ "$engine" == "blast" ]]; then
  if [ ! -f "$WD/blastdb/$GNM_TO_BASE".nin ]; then
    cat_or_zcat "$WD/genome_to/$GNM_TO_BASE" | 
      makeblastdb -in - -dbtype nucl \
                  -title "$GNM_TO_BASE" \
                  -hash_index \
                  -out "$WD/blastdb/$GNM_TO_BASE"
  fi
  
  if [ ! -f "$WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE".bln ]; then
    echo
    echo "== Run BLAST"
    NOW=$(date)
    echo "TIME START BLAST: $NOW"
    cat "$WD/marker_from/$MRK_FR_BARE".1kflank.fna | 
      blastn -db "$WD/blastdb/$GNM_TO_BASE" \
             -query - \
             -num_threads "$NPROC" \
             -evalue "$evalue" \
             -perc_identity "$perc_identity" \
             -outfmt "6 std qlen qcovs" |
               cat > "$WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bln"
    echo "DONE with BLAST"
    NOW=$(date)
    echo "TIME END BLAST:   $NOW"
  else 
    echo "== Skipping BLAST, as the output file exists."
  fi
elif [[ "$engine" == "burst" ]]; then
  if [ ! -f "$WD/blastdb/$GNM_TO_BASE.edx" ]; then
      burst_linux_DB12 -r "$WD/genome_to/$GNM_TO_BASE" \
                       -a "$WD/blastdb/$GNM_TO_BASE.acx" \
                       -o "$WD/blastdb/$GNM_TO_BASE.edx" \
                       -d DNA 1000 -i "0.$perc_identity" -s
  fi
  
  if [ ! -f "$WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bst" ]; then
    echo
    echo "== Run BURST"
    NOW=$(date)
    echo "TIME START BURST: $NOW"
    echo "burst_linux_DB12 -q $WD/marker_from/$MRK_FR_BARE.1kflank.fna -a $WD/blastdb/$GNM_TO_BASE.acx \\"
    echo "                 -r $WD/blastdb/$GNM_TO_BASE.edx -o $WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bst --mode FORAGE --forwardreverse"

    burst_linux_DB12 -q "$WD/marker_from/$MRK_FR_BARE.1kflank.fna" \
                     -a "$WD/blastdb/$GNM_TO_BASE.acx" \
                     -r "$WD/blastdb/$GNM_TO_BASE.edx" \
                     --mode FORAGE \
                     --forwardreverse \
                     -o "$WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bst"
    echo "DONE with BURST"
    NOW=$(date)
    echo "TIME END BURST:   $NOW"
  else 
    echo "== Skipping BURST, as the output file exists."
  fi
else 
  echo "The search engine must be specified as either \"burst\" (default) or \"blast\"."
  echo "The value of \"engine\" is currently set as $engine. Please check and correct the config file."
fi

echo
echo "== Filter $engine output and write new marker file (as a tsv file)"

if [[ "$engine" == "blast" ]]; then
  echo "$WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bln"
  # The blastn used above reports 14 fields, ending with with "qlen qcovs"
  cat "$WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bln" | top_line.awk | 
    sort -k2,2 -k1r,1r > "$WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bln.sorted"

  cat "$WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bln.sorted" |
    marker_blast_to_gff.pl -genome "$WD/genome_to/$GNM_TO_BASE" \
                           -gff_source "$gff_source" \
                           -gff_type "$gff_type" \
                           -gff_ID_prefix "$gff_ID_prefix" \
                           -max_len "$max_len" \
                           -qcov_identity "$qcov_identity" \
                           -sample_len "$sample_len" \
                           -gff_prefix_regex "$gff_prefix_regex" \
                           -out "$WD/marker_to/$marker_to" \
                           -verbose
elif [[ "$engine" == "burst" ]]; then
  # Calculate and add query sequence length and qcovs (percentage of the query that matches the target)

  # Note that in burst tabular output, the starting coordinate of the match is one less than would be reported by blast.
  # Also, the query length and the percent of the alignment relative to the query need to be added manually (cols 13 and 14).

  echo "$WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bst"
  seqlen "$WD/marker_from/$MRK_FR_BARE.1kflank.fna" | sort -k1,1 | uniq > "$WD/marker_from/$MRK_FR_BARE.len"
  join <(sort -k1,1 -k11n,11n "$WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bst") "$WD/marker_from/$MRK_FR_BARE.len" |
    perl -pe 's/ +/\t/g' | 
    awk '{print $0 "\t" int(100*($13-$11))/$13}' | 
    awk -v OFS="\t" '{ $9 = $9 + 1; print }' | sort -k1,1 |
    top_line.awk | sort -k2,2 -k1r,1r > "$WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bst.sorted"

  cat "$WD/blastout/$MRK_FR_BARE.x.$GNM_TO_BASE.bst.sorted" |
    marker_blast_to_gff.pl -genome "$WD/genome_to/$GNM_TO_BASE" \
                           -gff_source "$gff_source" \
                           -gff_type "$gff_type" \
                           -gff_ID_prefix "$gff_ID_prefix" \
                           -max_len "$max_len" \
                           -qcov_identity "$qcov_identity" \
                           -sample_len "$sample_len" \
                           -gff_prefix_regex "$gff_prefix_regex" \
                           -out "$WD/marker_to/$marker_to" -verbose
fi

echo
echo "== Sort GFF"
sort_gff.pl "$WD/marker_to/$marker_to.gff3" > "$WD/marker_to/tmp.gff"
mv "$WD/marker_to/tmp.gff" "$WD/marker_to/$marker_to.gff3"

echo
echo "== Compare the initial and mapped markers and report"
echo "==   Marker list 1: $WD/$work_dir/marker_from/lis.$MRK_FR_BARE"
echo "==   Marker list 2: $WD/$work_dir/marker_to/lis.$marker_to"
echo "==   Marker report: $WD/$work_dir/marker_to/report.${MRK_FR_BARE}--${marker_to}.tsv"

# Extract ID and allele from the "from" bed file, and print "+" orientation for all
cat "$WD/marker_from/$MRK_FR_BARE.bed" | awk -v OFS="\t" '{print $4, "+", $5}' | sort > "$WD/marker_from/lis.$MRK_FR_BARE"

# Next: Extract ID, allele, and orientation from the "to" bed file
cut -f4,6,7 "$WD/marker_to/$marker_to.bed" | sort > "$WD/marker_to/lis.$marker_to"

# Fields in joined result: markerID, orient, allele_gnm1, orient, allele_gnm2
#                              0        1  2  3  4
#                          ss715578401  +  G  +  G
#                          ss715578490  +  G  -  C
join -a1 "$WD/marker_from/lis.$MRK_FR_BARE" "$WD/marker_to/lis.$marker_to" |
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
echo "==   $WD/$work_dir/marker_to/$marker_to.bed"
echo "==   $WD/$work_dir/marker_to/$marker_to.gff3"
echo
echo "== Run completed. Look for results at $work_dir/marker_to/"
echo

exit 0

