#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Bio::Seq;
use Bio::DB::Fasta;
use feature "say";

my $usage = <<EOS;
  Synopsis: cat BLAST-OUTPUT | top_line.awk | marker_blast_to_gff.pl -genome GENOME_FILE -out FILE_BASE [options]

   filter_marker_blast_data.awk - Given input data consisting of blastn in tabular -m8 format,
   in which the query IDs have a suffix of ".UP" or ".DN", e.g.,
     glyma.Wm82.gnm1.ss715579418.UP
     glyma.Wm82.gnm1.ss715579418.DN
   report the position (molecule and basepair) immediately prior to the match
   of the ".DN" sequence.
  
   Note 1: The tabular BLAST output should be augmented with two additional fields, "qlen qcovhsp":
     -outfmt "6 std qlen qcovhsp"

   Note 2: The .UP and .DN matches represent upstream and downstream flanking sequences
     adjacent to a marker; so the position prior to the .DN match should be the
     marker location in the target sequence if the .UP and .DN sequences matched the
     target in forward orientation; and the position prior to the .UP match should be the
     marker location in the target sequence if the .UP and .DN sequences matched the
     target in reverse orientation.

   The BLAST output should be filtered such that the matches of the flanking query sequences
   (.UP and .DN) are near-perfect top matches, and sorted relative to the target sequence and
   such that the .UP and .DN for a given marker are adjacent in the filtered BLAST output.
   The standard ordering of the BLAST output should be sufficient.

  Required:
    STDIN with BLAST output, created with -outfmt "6 std qlen qcovhsp"
    -genome      Target genome file in which variants are to be mapped (uncompressed)
    -out         Name for new marker files, sans extension. Three files will be written: .bed, .gff3, .log
    
  Options:   
    -qcov_identity  Percent identity for query coverage (qcovhsp), in range 0..100 [90]
    -sample_len  Maximum length of sequence variant to report, as a sample, in the GFF 9th column [10]
    -max_len     Maximum variant length for which to report a GFF line [200]
                   For a SNP, most variants are 1 base long, but there may be longer indels. 
                   SSRs are expected to be longer, but very long matches probably indicate either
                   poor or repetitive matches of the flanking sequences, or large genomic rearrangements.
    -gff_source  String to use for the gff3 source (column 2); typically, a project, data source, or program
    -gff_type    String to use for the gff3 type (column 3), e.g. "genetic_marker" or other sequence ontology term
    -gff_prefix_regex  Regular expression for stripping the genome prefix from the marker IDs, e.g. glyma.Wm82.gnm1.Sat_413 => Sat_413
    -verbose     (boolean) Some debug output, marked with "=="
    -help        (boolean) This message.
EOS

my ($genome_file, $out_file, $verbose, $help);
my $qcov_identity=90;
my $sample_len=10;
my $max_len=200;
my $gff_source="map-markers-to-assembly";
my $gff_type="genetic_marker";
my $gff_ID_prefix="TARGET_GENOME_VERSION";
my $gff_prefix_regex="^[^.]+\.[^.]+\.[^.]+\.";

GetOptions (
  "genome=s"     => \$genome_file,   
  "out=s"        => \$out_file,   
  "qcov_identity:i"   => \$qcov_identity,
  "sample_len:i" => \$sample_len,
  "max_len:i"    => \$max_len,
  "gff_source:s" => \$gff_source,
  "gff_type:s"   => \$gff_type,
  "gff_ID_prefix:s" => \$gff_ID_prefix,
  "gff_prefix_regex:s" => \$gff_prefix_regex,
  "verbose"      => \$verbose,
  "help"         => \$help,
);

die "$usage" if $help;
die "$usage\nPlease specify a GFF output filename with -out\n" unless $out_file;
die "$usage\nPlease specify a genome assembly with -genome\n" unless $genome_file;
die "Can't open in $genome_file: $!\n" unless (-f $genome_file);

my $gff_file_base = $out_file;
$gff_file_base =~ s/\.gff3?$//;

my ($GFF_FH, $BED_FH, $LOG_FH);
open ( $GFF_FH, ">", "$gff_file_base.gff3" ) or die "Can't open out $gff_file_base.gff3: $!\n";
open ( $BED_FH, ">", "$gff_file_base.bed") or die "Can't open out $gff_file_base.bed: $!\n";
open ( $LOG_FH, ">", "$gff_file_base.log") or die "Can't open out $gff_file_base.log: $!\n";

my $REX;
if ($gff_prefix_regex){ $REX=qr/$gff_prefix_regex/ }
else { $REX=qr/$/ }

my ($prev_UD, $prev_base, $prev_start, $prev_end);
my %seen_markID;
my %seen_skippedID;

# Load genome into bioperl object
my $db = Bio::DB::Fasta->new($genome_file);

# Read in the BLAST output and process
while (<>) {
  chomp;
  next if /^#/;
  my @F = split /\t/;
  my ($markID_UD, $seqID, $this_start, $this_end, $bitsc, $qcovhsp) = 
    ($F[0], $F[1], $F[8], $F[9], $F[11], $F[13]);
  $markID_UD =~ /(.+)\.(UP|DN)$/;
  my ($base_id, $up_or_dn) = ($1, $2);

  # Strip prefix from ID from the "from" markers
  my $name=$base_id;
  $name =~ s/$REX//;

  unless ($seen_markID{$name}){ $seen_markID{$name}++; }

  if ($qcovhsp < $qcov_identity){
    if ($verbose){
      say "== Skipping $name $up_or_dn bbecause qcovhsp<qcov_identity: $qcovhsp<$qcov_identity";
    }
    say $LOG_FH "Skipping $name $up_or_dn because qcovhsp<qcov_identity: $qcovhsp<$qcov_identity";
    unless ($seen_skippedID{$name}){ $seen_skippedID{$name}++; }
  }
  else {
    if (defined $prev_base && $base_id eq $prev_base) {
      if ($this_start > $prev_end && $up_or_dn =~ /DN$/) { # forward-forward
        if ($verbose){say "AA: handling forward-forward; $this_start > $prev_end && $up_or_dn" }
        my ($short_var, $full_var);
        if ( ($this_start-1)-$prev_end < 0 ) {
          my ($start, $end) = ($prev_end, $this_start-1);
          say $LOG_FH "Skipping $name because start is greater than end: $start, $end";
          next;
        }
        elsif ( ($this_start-1)-$prev_end == 0 ) { # marker is of zero length - probably an indel
          ($short_var, $full_var) = ("_", "_");
          if ($verbose){say "AA: $this_start-1, $prev_end, $short_var";}
          $prev_end--; # adjust the start coord to preserve the (collapsed) marker
        }
        else { # start-end (bed coords) are >=1
          ($short_var, $full_var) = get_variant($seqID, $prev_end + 1, $this_start - 1, "FWD");
        }
        if ($short_var =~ /WARN/){
          if ($verbose){
            say "== Skipping $name $short_var";
          }
          say $LOG_FH "Skipping $name $short_var";
          unless ($seen_skippedID{$name}){ $seen_skippedID{$name}++; }
          next;
        }
        else {
          my $ninth = "ID=$gff_ID_prefix$name;Name=$name;ref_allele=$short_var";
          unless ($seen_skippedID{$name}){
            say $GFF_FH join("\t", $seqID, $gff_source, $gff_type, $prev_end + 1, $this_start - 1, ".", ".", "+", $ninth );
            say $BED_FH join("\t", $seqID,                         $prev_end,     $this_start - 1, $name, $bitsc, "+", $full_var);
          }
        }
      } 
      elsif ($this_start <= $prev_end && $up_or_dn =~ /DN$/) { # forward-reverse
        if ($verbose){say "BB: handling forward-reverse; $this_start <= $prev_end && $up_or_dn"}
        my ($short_var, $full_var);
        if ( ($prev_end-1)-$this_start < 0 ) {
          my ($start, $end) = ($this_start, $prev_end-1);
          say $LOG_FH "Skipping $name because start is greater than end: $start, $end";
          next;
        }
        elsif ( ($prev_end-1)-$this_start == 0 ) { # marker is of zero length - probably an indel
          ($short_var, $full_var) = ("_", "_");
          if ($verbose){say "BB: $this_start-1, $prev_end, $short_var";}
          $this_start--; # adjust the start coord to preserve the (collapsed) marker
        }
        else { # start-end (bed coords) are >=1
          ($short_var, $full_var) = get_variant($seqID, $this_start + 1, $prev_end - 1, "REV"); # seq here is reverse-complemented
        }
        if ($short_var =~ /WARN/){
          if ($verbose){
            say "== Skipping $name $short_var";
          }
          say $LOG_FH "Skipping $name $short_var";
          unless ($seen_skippedID{$name}){ $seen_skippedID{$name}++; }
          next;
        }
        else {
          my $ninth = "ID=$gff_ID_prefix$name;Name=$name;ref_allele=$short_var";
          unless ($seen_skippedID{$name}){
            say $GFF_FH join("\t", $seqID, $gff_source, $gff_type, $this_start + 1, $prev_end - 1, ".", ".", "-", $ninth );
            say $BED_FH join("\t", $seqID,                         $this_start,     $prev_end - 1, $name, $bitsc, "-", $full_var);
          }
        }
      }
    }
  }

  if (!defined $prev_base || $base_id ne $prev_base || $. == 1) {
    $prev_UD = $up_or_dn;
    $prev_base = $base_id;
    $prev_start = $F[8];
    $prev_end = $F[9];
  }
}

my $ct_markers = scalar(keys %seen_markID);
my $ct_skipped = scalar(keys %seen_skippedID);
say "== Of $ct_markers markers, $ct_skipped were skipped as not meeting specified parameters.";
say "== See log file at $gff_file_base.log";

#####################

sub get_variant {
  my ($id, $start, $end, $orient) = @_;
  my $full_seq;
  if ($orient =~ /FWD/){
    $full_seq = $db->seq($id, $start, $end);
  }
  else { # $orient is REV
    $full_seq = $db->seq($id, $start, $end);
    my $seqobj = Bio::Seq->new(-seq => $full_seq, -alphabet => "dna");
    my $revcom_obj = $seqobj->revcom();
    my $revcom_seq = $revcom_obj->seq();
    #say "SEQ: " . substr($full_seq, 0, $max_len);
    #say "REV: " . substr($revcom_seq, 0, $max_len);
    $full_seq = $revcom_seq;
  }
  my $seq = $full_seq;
  my $len = length($full_seq);
  if ( $len > $sample_len && $len <= $max_len ){
    $seq = substr($full_seq, 0, $sample_len);
    $seq = $seq . "_from_$len" . "_nt";
  }
  if ( $len > $sample_len && $len > $max_len ){
    $seq = substr($full_seq, 0, $max_len);
    $seq = "WARN variant is $len" . " nt. First $max_len nt: " . $seq;
  }
  return($seq, $full_seq);
}
__END__

Steven Cannon
2025-01-29 First functional version.
2025-01-31 Require GFF output filename and also print to a derived bed file.
2025-02-04 Print warnings to a log file.
2025-02-05 Add orientation and score to BED output, and report rev-complimented sequence if mapping to negative strand
2026-01-14 Change handling of identity, now calling it qcov_identity to distinguish it from percent identity
