#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max);
use Bio::Seq;       
use Bio::DB::Fasta; 
use feature "say";

my $usage = <<EOS;
  Synopsis: marker_gff_to_bed_and_var.pl -markers MARKER_FILE -genome GENOME_FILE -out OUT_FILE [options]
 
  Read marker coordinates in GFF format (compressed or uncompressed) and 
  a genome file (uncompressed), with samtools faidx index.
  Return a four-column BED file with genomic coordinates flanking the marker coordinates,
  plus the sequence (presumably the marker variant) under those coordinates.
  Adjust from GFF 1-based, closed [start, end] to BED 0-based, half-open [start-1, end)
  Unless the "-merge" flag is set, two sequences are returned for each marker: 
  one before (>ID.UP), one after (>ID.DN). 

  The fai file (or a file with seqID in column 1 and molecule length in column 2) 
  is used to identify markers near the ends of molecules.
  The index file should be in the same directory as the genome file, as FILENAME.fai.
  The optional -pad flag indicates the distance up- and down-stream from the marker
  to be extracted (default 1000). If the marker position is less than (pad) units from 
  the start of a molecule, then the padded start will be set to 0; and if the marker
  position is less than (pad) units from the end of the molecule, then the padded end
  will be set to the length of that molecule.
  
  Required:
    -markers GFF file (either compressed or not)
    -genome  Genome file corresponding with GFF file. Should be uncompressed, and with accompanying index file.
    -out     Name for output file (bed), sans extension. Two files will be written: FILE.bed and FILE.UD.bed
               FILE.bed    with  seqID, marker name, coordinates, and allele sequence
               FILE.UD.bed with  seqID, marker name.UP, coordinates
                                 seqID, marker name.DN, coordinates
       
  Options:
    -pad   The flanking distance up- and down-stream from the marker. Default 1000.
    -min_flank   Minimum length of UP or DN sequence to report. Default 100.
    -gff_prefix_regex  Regular expression for stripping the genome prefix from the marker IDs, e.g. glyma.Wm82.gnm1.Sat_413 => Sat_413
    -merge (boolean; default false). If set, then one sequence will be returned for each marker;
             otherwise, return one UP (>ID.UP), one DNer (>ID.DN).
    -verbose (boolean) Some debug output, marked with "=="
    -help  (boolean) This message.
EOS

my ($genome_file, $marker_file, $out_file, $merge, $verbose, $help);
my $pad = 1000;
my $min_flank = 100;
my $gff_prefix_regex="^[^.]+\.[^.]+\.[^.]+\.";

GetOptions (
  "genome=s"    =>  \$genome_file,   
  "markers=s"   =>  \$marker_file,   
  "out=s"       =>  \$out_file,   
  "pad:i"       =>  \$pad,  
  "merge"       =>  \$merge,
  "gff_prefix_regex:s" => \$gff_prefix_regex,
  "min_flank:i" => \$min_flank,
  "verbose"     =>  \$verbose,
  "help"        =>  \$help,
);

die "$usage\nPlease specify a genome assembly with -genome\n" unless $genome_file;
die "$usage\nPlease specify a marker file with -marker\n" unless $genome_file;
die "$usage\nPlease specify an output filename (sans extension) with -out\n" unless -e "$genome_file.fai";
die "$usage" if $help;
die "Can't open in $genome_file: $!\n" unless (-f $genome_file);

open my $FAI_FH, "<", "$genome_file.fai" or die "Can't open in $genome_file.fai: $!\n";

# Strip extension in case one was specified (against instructions)
my $mrk_file_base = $out_file;
$mrk_file_base =~ s/\.bed$//;
$mrk_file_base =~ s/\.gff3?$//;

my ($BED_FH, $BED_UD_FH);
open ( $BED_FH, ">", "$mrk_file_base.bed" ) or die "Can't open out $mrk_file_base.bed$!\n";
open ( $BED_UD_FH, ">", "$mrk_file_base.UD.bed") or die "Can't open out $mrk_file_base.UD.bed: $!\n";

my $REX;
if ($gff_prefix_regex){ $REX=qr/$gff_prefix_regex/ }
else { $REX=qr/$/ }

# Read in fai file and get molecule lengths
my %seqID_len;
while (<$FAI_FH>) {
  chomp;
  my @fields = split(/\t/, $_);
  my ($seqID, $len) = ($fields[0], $fields[1]);
  $seqID_len{$seqID} = $len;
  #if ($verbose){say "== AA: $seqID\t$len";}
}

# Load genome into bioperl object
my $db = Bio::DB::Fasta->new($genome_file);

# Read in the sequence
my $MRK_FH;
if ( $marker_file =~ /gz$/ ){
  open( $MRK_FH, "zcat $marker_file|" ) or die "Can't do zcat $marker_file| : $!";
}
else {
  open ( $MRK_FH, "<", $marker_file ) or die "Can't open in $marker_file: $!\n";
}

# Read in the GFF;
while (<$MRK_FH>) {
  s/\r?\n\z//; # CRLF to LF
  chomp;
  my $line = $_;
  if ( $line =~ /^#|^\s*$/ ) { next } # Skip comment and spacer lines
  else { # body of the GFF
    my @parts = split(/\t/, $line);
    if (scalar(@parts)<9){ next }
    my $col9 = $parts[8];
    my @col9_attrs = split(/;/, $col9);
    my ($seqID, $mrk_start, $mrk_end) = ($parts[0], $parts[3]-1, $parts[4]); # in BED coords
    my $mrk_id;
    for my $attr (@col9_attrs){
      if ($attr =~ /ID=(.+)/){
        $mrk_id = $1;
        $mrk_id =~ s/([^;]+);.+/$1/;
      }
    }
    my $pad_start = max(0, $mrk_start-$pad);

    if ($verbose){
      say "== BB: mrk_id:\t$mrk_id";
      say "== BB: seqID:\t$seqID";
      say "== BB: mrk_sta:\t$mrk_start";
      say "== BB: mrk_end:\t$mrk_end";
      say "== BB: pad:\t$pad";
      say "";
    }

    my $pad_end = min($mrk_end+$pad, $seqID_len{$seqID});
    my $variant = $db->seq($seqID, $mrk_start+1, $mrk_end);

    # Strip prefix from ID from the marker IDs in the bed file with variants
    my $name=$mrk_id;
    $name =~ s/$REX//;

    my $size_left = $mrk_start-$pad_start;
    my $size_right = $pad_end-$mrk_end;
    if ( $mrk_start-$pad_start<$min_flank || $pad_end-$mrk_end<$min_flank ) {
      if ($verbose){
        say "== CC: Skipping $seqID because flanking seq is too short: $mrk_start-$pad_start or $pad_end-$mrk_end";
        say "== CC: $mrk_start-$pad_start, $pad_end-$mrk_end, $min_flank";
        say "== CC: $size_left, $size_right, $min_flank";
        say "";
      }
      next;
    }
    else {
      if ($merge){ # Return one sequence per marker, including flanking upstream and down (in BED coords)
        say $BED_UD_FH join("\t", $seqID, $pad_start, $pad_end, $mrk_id);

        say $BED_FH    join("\t", $seqID, $pad_start, $pad_end, $name, $variant);
      }
      else { # Return two sequences per marker: one upstream, one down (in BED coords)
        say $BED_UD_FH join("\t", $seqID, $pad_start, $mrk_start, "$mrk_id.UP");
        say $BED_UD_FH join("\t", $seqID, $mrk_end, $pad_end, "$mrk_id.DN");

        say $BED_FH    join("\t", $seqID, $mrk_end, $pad_end, $name, $variant);
      }
    }
  }
}

__END__

Steven Cannon
2023-03-09 Start 
2023-11-28 Add -verbose flag and some debug reporting
2015-01-22 Remove 1-padding from $mrk_end. Add option for reporting variant sequence
2025-01-28 Make reporting of variant sequence non-optional
2025-01-31 Report two types of BED files -- one with .UP and .DN for each marker, and one with the seq variant
2026-01-13 Add parameter min_flank and suppress printing sequences that are below this value
