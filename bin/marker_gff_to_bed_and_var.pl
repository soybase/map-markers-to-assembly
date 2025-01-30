#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max);
use Bio::Seq;       
use Bio::DB::Fasta; 
use feature "say";

my $usage = <<EOS;
  Synopsis: cat GFF_FILE.gff3 | marker_gff_to_bed_and_var.pl -genome GENOME_FILE [options]
  
  Read marker coordinates in GFF format and a genome file (uncompressed), with samtools faidx index.
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
    stream of GFF data on STDIN
    -genome  Genome file corresponding with GFF file. Should be uncompressed, and with accompanying index file.
       
  Options:
    -pad   The flanking distance up- and down-stream from the marker. Default 1000.
    -merge (boolean; default false). If set, then one sequence will be returned for each marker;
             otherwise, return one beUP (>ID.UP), one DNer (>ID.DN).
    -out   File to write to; otherwise, to stdout.
    -verbose (boolean) Some debug output, marked with "=="
    -help  (boolean) This message.
EOS

my ($genome_file, $out_file, $merge, $verbose, $help);
my $pad = 1000;

GetOptions (
  "genome=s" =>  \$genome_file,   
  "pad:i"    =>  \$pad,  
  "merge"    =>  \$merge,
  "out:s"    =>  \$out_file,   
  "verbose"  =>  \$verbose,
  "help"     =>  \$help,
);

die "$usage" unless $genome_file;
die "$usage" unless -e "$genome_file.fai";
die "$usage" if $help;
if ($genome_file){
  die "Can't open in $genome_file: $!\n" unless (-f $genome_file);
}

open my $FAI_FH, "<", "$genome_file.fai" or die "Can't open in $genome_file.fai: $!\n";

my $OUT_FH;
if ( $out_file ){
  open ( $OUT_FH, ">", $out_file ) or die "Can't open out $out_file: $!\n";
}
else {
  open ( $OUT_FH, ">&", \*STDOUT) or die;
}

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

# Read in the GFF;
while (<>) {
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
    if ($merge){ # Return one sequence per marker, including flanking upstream and down (in BED coords)
      say $OUT_FH join("\t", $seqID, $pad_start, $pad_end, $mrk_id);
    }
    else { # Return two sequences per marker: one upstream, one down (in BED coords)
      my $variant = $db->seq($seqID, $mrk_start+1, $mrk_end);
      say $OUT_FH join("\t", $seqID, $pad_start, $mrk_start, "$mrk_id.UP", $variant);
      say $OUT_FH join("\t", $seqID, $mrk_end, $pad_end, "$mrk_id.DN", $variant);
    }
  }
}

__END__

Steven Cannon
2023-03-09 Start 
2023-11-28 Add -verbose flag and some debug reporting
2015-01-22 Remove 1-padding from $mrk_end. Add option for reporting variant sequence
2025-01-28 Make reporting of variant sequence non-optional
