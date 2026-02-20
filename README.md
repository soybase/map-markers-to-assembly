# map-markers-to-assembly

## Intro and usage
This little workflow is driven by the shell script map-markers.sh.

It extracts flanking sequence around a sequenced-based marker from one genome assembly, then searches 
for the best corresponding sequence in a second assembly, and reports the locations of the SNP 
marker in the second assembly.

The output, which goes to the "work_dir/marker_to" directory, includes three report files:

  * orient.FROM_GNM--TO_GNM.txt
  * report.FROM_GNM--TO_GNM.tsv
  * TO_FILE.log

... and three primary outputs:
  * TO_GNM_abs.bed   # in molecule (absolute/FWD) orientation
  * TO_GNM_rel.bed   # in relative orientation (relative between the FROM and TO genomes)
  * TO_GNM.gff3

Note that the markers and alleles are all reported in the molecule (FWD) orientation in the GFF
and in the `_abs.bed` file, whereas in the `_rel.bed` file, the alleles are reported for the 
"to" genome relative to the "from" genome around those markers. In the "report" file, the alleles 
are given both in the target FWD orientation and in the relative orientation:
```
                                              from
                                                     to;BED
                                                            to;GFF
                                              vvv    vvv    vvv
  #markerID  compare   len1   len2   orient   var1   var2   var2_fwd
  snp1       same      1      1      +        T      T      T
  snp2       same      1      1      +        G      G      G
  snp3       NOT       1      1      +        G      A      A
  snp4       same      1      1      -        T      T      A
  snp5       same      1      1      -        T      T      A
  snp6       same      1      1      -        G      G      C
  snp7       same      1      1      -        A      A      T
  snp8       NOT       1      2      -        C      CT     AG
  snp9       NOT       1      1      -        C      T      A
```

The following are the main steps:

  * Create an index on the uncompressed FROM genome;
  * Put the marker information into four-column BED format, with 1000 bases on each side of the SNP;
  * Extract the sequences from the FROM genome;
  * Run BLAST or burst (depending on the chosen search engine) against the TO genome;
  * Filter BLAST or burst output and writes new marker file (as a gff3 file).

The dependencies (see installation instructions below) are:
  * samtools 
  * bedtools 
  * BioPerl
  * blast or burst (blast preferred)

and eight scripts in the bin directory:
  * map-markers.sh
  * marker_blast_to_gff.pl
  * marker_gff_to_bed_and_var.pl
  * marker_report.pl
  * top_line.awk
  * sort_gff.pl
  * filter_fasta_for_Ns.awk
  * liftover_vcf.pl

```bash
NAME
  map-markers-to-assembly.sh  -- Given a file of marker locations in one genome, report the
    locations of those markers in a second genome.

SYNOPSIS
  map-markers.sh  -c CONFIG_FILE

  Required:
           -c (path to the config file)

  Options: -h help

```

The script `liftover_vcf.pl` is a utility that can be run after `map-markers.sh` to update coordinates and
alleles in a VCF file, projecting from the `$from_genome` to the `$to_genome`. The script use as an input 
the `${marker_to}_abs.bed` file, like so:
```
    ./bin/liftover_vcf.pl --vcf variants.vcf.gz \
                          --bed work_dir/marker_to/${marker_to}_abs.bed \
                          --reference $to_genome \
                          --output lifted_variants.vcf
```
The `liftover_vcf.pl` script is quick, so can be run in interactive mode rather than via SLURM batch script; 
however, if running it interactively, you will need to activate the `map-markers` conda environment (see details below).
Or it could be run via a batch script.

## Config file
Specify paths to the "from" and "to" genome assemblies and the (gff3) marker file.
The two genome files will be uncompressed in the work directory, as part of this script. The primary intended
use of the script is on the same file system as the data files, in which case the starting files can be copied
into the work directory. If they are coming from a remote remote location, one solution would be to pull the
files locally into a data directory using scp or equivalent, and then give the paths to those files.

  VARIABLES set in config file:
``` bash
    marker_from   - Full filepath to file with marker names and locations on first Genome; in gff3 format, compressed
    genome_from   - Full filepath to first genome assembly, corresponding with the coordinates in the marker_from file; compressed
    genome_to     - Full filepath to the second genome assembly, to which the markers will be projected
    marker_to     - Name for new marker files, sans extension. Three files will be written: FILE.bed, FILE.gff3, FILE.log
    gff_source    - String to use for the gff3 source (column 2); typically, a project, data source, or program
    gff_ID_prefix - String to use for the ID prefix in the gff3 9th column -- for example, the genome assembly version
```

... and some config-file parameters that can generally be left with default values:
```
    engine        - blast or burst [blast]
                    BLAST should work well in essentially every situation. The reason to consider BURST is that it is
                    much faster for very large (100k+) marker sets. The downsides to BURST are lower sensitivity
                    (~5% vs. BLAST in this context), and the target-genome index files are about 20x larger than
                    for BLAST; and their creation takes a large amount of memory (500 GB for a typical genome).

    gff_type      - String to use for the gff3 type (column 3), e.g. "genetic_marker" or other SOFA sequence ontology term
    gff_prefix_regex - Regular expression for stripping the genome prefix from the FROM marker IDs, e.g. glyma.Wm82.gnm1.Sat_413 => Sat_413
    perc_identity - Minimum percent identity in range 0..100 for BLAST or BURST percent identity [95]
    qcov_identity - Minimum percent identity in range 0..100 for BLAST or BURST qcovhsp [80]
    sample_len    - Maximum length of sequence variant to report, as a sample, in the GFF 9th column [10]
    max_var_len   - Maximum variant length for which to report a GFF line [25]
    min_flank     - Minimum length of flanking around a marker, to be used in search of the target genome. [100]
    work_dir      - Work directory; [work_dir]
```

In a HPC environment, the workflow should be called using a job manager -- either interactively or via a job submission script.
Mapping a few hundred markers should take only a few minutes, but mapping tens or hundreds of thousands of markers may take 
on the order of an hour.

### Installation of dependencies:

Create a conda environment (unless an environment with the necessary packages exists). 
The conda environment specified by the `environment.yml` config is called `map-markers`.

    `conda env create`

NOTE for the SoyBase and Legume Information System project: a suitable project onvironment is `ds-curate`; in that case, just do `source activate ds-curate`.
See details regarding the [ds-curate environment here](https://github.com/legumeinfo/datastore-specifications/tree/main/PROTOCOLS#set-paths-and-start-a-conda-environment-with-software-needed-for-curation).

## Running the program:

If running the program interactively, start an interactive session, e.g. 
```
  salloc  
  module load miniconda
  source activate map-markers
  PATH=$PWD/bin:$PATH

  map-markers.sh -c config/gnm1_to_gnm2_SoySSR.conf
```

Or if running the program through a job control script, modify `batch_map_markers_example.sh`, which 
does essentially the same things as in the interactive invocation, but via script.

