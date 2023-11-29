#!/usr/bin/awk -f
#
# NAME
#   filter_marker_blast_data.awk - Given input data consisting of blastn in tabular -m8 format,
#     in which the query IDs have a suffix of ".UP" or ".DN", e.g.,
#       glyma.Wm82.gnm1.ss715579418.UP
#       glyma.Wm82.gnm1.ss715579418.DN
#     report the position (molecule and basepair) immediately prior to the match 
#     of the ".DN" sequence.
#
#     The .UP and .DN matches represent upstream and downstream flanking sequences
#     adjacent to a marker; so the position prior to the .DN match should be the
#     marker location in the target sequence if the .UP and .DN sequences matched the 
#     target in forward orientation; and the position prior to the .UP match should be the
#     marker location in the target sequence if the .UP and .DN sequences matched the 
#     target in reverse orientation.
# 
#     The BLAST output should be filtered such that the matches of the flanking query sequences
#     (.UP and .DN) are near-perfect top matches, and sorted relative to the target sequence and
#     such the .UP and .DN for a given marker are adjacent in the filtered BLAST output. 
#
# SYNOPSIS
#   cat BLAST-OUTPUT | awk '$4>=90 && $14>=90' | top_line.awk | sort -k2,2 -k9n,9n | filter_marker_blast_data.awk
#     In this filtering pattern, the standard BLAST output is augmented with three additional fields:
#       -outfmt "6 std qlen qcovs qcovhsp"
#     so $4>=90 requires >=90% identity and $14>=90 requires >=90% query-coverage-per-sequence.
#     
# AUTHOR  Steven Cannon

BEGIN { FS = OFS = "\t" }

base_id = substr($1,1,length($1)-3) { }
up_or_dn = substr($1,length($1)-1,2) { }
base_id == prev_base && prev_UD eq "UP" && $1 ~ /DN$/ { print base_id, $2, $9-1, "FWD"; prev_UD=up_or_dn; prev_base=base_id }
base_id == prev_base && prev_UD eq "DN" && $1 ~ /UP$/ { print base_id, $2, $10-1, "REV"; prev_UD=up_or_dn; prev_base=base_id }
base_id != prev_base || NR==1 { prev_UD=up_or_dn; prev_base=base_id }

