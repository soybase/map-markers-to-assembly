# map-markers-to-assembly

This little workflow is driven by the shell script map-markers.sh.

It extracts flanking sequence around a SNP marker from one genome assembly, then searches for the 
best corresponding sequence in a second assembly, and reports the locations of the SNP 
marker in the second assembly.
Additionally, the process reports the reference allele for each variant in the FROM and the TO genomes.

The following are the main steps:

  * Create an index on the uncompressed FROM genome;
  * Put the marker information into four-column BED format, with 1000 bases on each side of the SNP;
  * Extract the sequences from the FROM genome;
  * Run BLAST against the TO genome;
  * Filter BLAST output and writes new marker file (as a gff3 file).

The dependencies are:
  * samtools 
  * bedtools 
  * BioPerl
  * blast+

and four scripts in the bin directory:
  * map-markers.sh
  * marker_blast_to_gff.pl
  * marker_gff_to_bed_and_var.pl
  * top_line.awk

```bash
NAME
  map-markers-to-assembly.sh  -- Given a file of marker locations in one genome, report the
    locations of those markers in a second genome.

SYNOPSIS
  map-markers.sh  -c CONFIG_FILE

  Required:
           -c (path to the config file)

  Options: -h help

VARIABLES set in config file:
  All input files should be compressed (bgzip). The two genome files will be uncompressed in the work directory,
  as part of this script. The primary intended use of the script is on the same file system as the data files,
  in which case the starting files can be copied into the work directory. If they are coming from a remote
  remote location, one solution would be to pull the files locally into a data directory using scp or equivalent.

    marker_from - Full filepath to file with marker names and locations on first Genome; in gff3 format, compressed

    genome_from - Full filepath to first genome assembly, corresponding with the coordinates in the marker_from file; compressed
    genome_to   - Full filepath to the second genome assembly, to which the markers will be projected

    marker_work - Name for intermediate/working marker file (bed), based on marker_from, sans extension. 
                    Two files will be written: FILE.bed and FILE.UD.bed

    marker_to   - Name for new marker files, sans extension. Two files will be written: FILE.bed and FILE.gff3

    gff_source  - String to use for the gff3 source (column 2); typically, a project, data source, or program
    gff_type    - String to use for the gff3 type (column 3), e.g. "genetic_marker" or other SOFA sequence ontology term
    gff_ID_prefix - String to use for the gff3 type (column 3), e.g. "genetic_marker" or other SOFA sequence ontology term

    identity    - Minimum percent identity in range 0..100 for blastn qcovhsp [90]
    sample_len  - Maximum length of sequence variant to report, as a sample, in the GFF 9th column [10]
    max_len     - Maximum variant length for which to report a GFF line [200]

    work_dir    - work directory; default work_dir
```

The workflow should generally not be especially computationally intensive, but in a HPC environment,
it should be called using a job manager -- either interactively or via a job submission script.

For example, interactively:

```bash
  salloc

  module load miniconda

  source activate ds-curate

  PATH=$PWD/bin:$PATH

  map-markers.sh -c config/gnm1_to_gnm2.conf 
```


