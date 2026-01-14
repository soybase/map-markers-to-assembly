#!/usr/bin/awk -f
#
# NAME
#   filter_fasta_for_Ns.awk - Suppress fasta sequences that have fewer than specified non-N bases.
#
# SYNOPSIS
#   cat FILE.fasta  | -v min_nonN=100 filter_fasta_for_Ns.awk
#
# AUTHOR: Steven Cannon

#!/usr/bin/env awk -f

BEGIN {
    # Default minimum number of non-N bases
    if (min_nonN == "")
        min_nonN = 100
}

/^>/ {
    # Process previous sequence
    if (seq != "") {
        nonN = seq
        gsub(/[Nn]/, "", nonN)
        if (length(nonN) >= min_nonN) {
            print header
            print seq
        }
    }
    header = $0
    seq = ""
    next
}

{
    seq = seq $0
}

END {
    # Process last sequence
    if (seq != "") {
        nonN = seq
        gsub(/[Nn]/, "", nonN)
        if (length(nonN) >= min_nonN) {
            print header
            print seq
        }
    }
}

