# map-markers-to-assembly

This little workflow is driven by the shell script map-markers-to-assembly.sh.

It extracts flanking sequence around a SNP marker from one genome assembly, then searches for the 
best corresponding sequence in a second assembly, and reports the locations of the SNP 
marker in the second assembly.

The following are the main steps:

  * Create an index on the uncompressed FROM genome;
  * Put the marker information into four-column BED format, with 1000 bases on each side of the SNP;
  * Extract the sequences from the FROM genome;
  * Run BLAST against the TO genome;
  * Filter BLAST output and writes new marker file (as a tsv file).

The dependencies are:
  * samtools 
  * bedtools 
  * blast+

and three scripts in the bin directory:
  * filter_marker_blast_data.awk
  * marker_gff_to_bed.pl
  * top_line.awk


OPERANDS
```
  Paths to four files (including the new marker-locations file to be created),
  and number of threads to use in blast search
    marker_locs    - Path to file of marker names and locations on Genome assembly 1; in gff3 format
    genome_from    - Path to genome assembly 1, corresponding with the coordinates in the marker_locs file
    genome_to      - Path to genome assembly 2, to which the markers will be projected
    marker_to_file - Path for the marker-locations file to be created; may include a directory path
    threads        - Number of threads to use in blast search
```

Note that the genome files will be uncompressed (gunzip) if they are in a compressed state,
as the bedtools getfasta command sometimes fails on compressed data.

In a multiprocessor machine with job management, the script should be called with 
a batch script, for example like the following:
```
  module load samtools bedtools blast+

  marker_locs=markers/FILE.gff3.gz
  genome_from=genomes/GENOME1/FILE.fna.gz
  genome_to=genomes/GENOME2/FILE.fna.gz
  marker_to_file=markers/FILE_NEW.tsv
  threads=20

  ./map-markers-to-assembly.sh $marker_locs $genome_from $genome_to $marker_to_file $threads
```
