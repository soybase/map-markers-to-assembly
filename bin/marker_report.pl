#!/usr/bin/env perl

use strict;
use warnings;

# Check command line arguments
if (@ARGV != 2) {
    print_usage();
    exit 1;
}

my ($file1, $file2) = @ARGV;

# Check if files exist
die "Error: Cannot find file '$file1'\n" unless -f $file1;
die "Error: Cannot find file '$file2'\n" unless -f $file2;

# Read and process the files
my %file1_data = read_file($file1);
my %file2_data = read_file($file2);

# Print header with new eighth field
print join("\t", "#markerID", "compare", "len1", "len2", "orient", "var1", "var2", "var2_fwd") . "\n";

# Process all markers from file1 (equivalent to join -a1)
foreach my $marker_id (sort keys %file1_data) {
    my $record1 = $file1_data{$marker_id};
    
    if (!exists $file2_data{$marker_id}) {
        # Marker not in target genome (equivalent to scalar(@F)==3)
        print join("\t", 
            $marker_id, 
            "NULL", 
            length($record1->{sequence}), 
            0, 
            ".", 
            $record1->{sequence}, 
            "NULL",
            "NULL"
        ) . "\n";
    } else {
        # Marker exists in both files
        my $record2 = $file2_data{$marker_id};
        my $orient = $record2->{orient};  # Use orientation from file2 (second field)
        my $var2_sequence = $record2->{sequence};
        
        # Calculate var2_fwd: reverse complement if orientation is negative
        my $var2_fwd;
        if ($orient eq '-') {
            $var2_fwd = reverse_complement($var2_sequence);
        } else {
            $var2_fwd = $var2_sequence;
        }
        
        if ($record1->{sequence} eq $record2->{sequence}) {
            # Allele is the same in both genomes
            print join("\t", 
                $marker_id, 
                "same", 
                length($record1->{sequence}), 
                length($record2->{sequence}), 
                $orient, 
                $record1->{sequence}, 
                $record2->{sequence},
                $var2_fwd
            ) . "\n";
        } else {
            # Alleles are different
            print join("\t", 
                $marker_id, 
                "NOT", 
                length($record1->{sequence}), 
                length($record2->{sequence}), 
                $orient, 
                $record1->{sequence}, 
                $record2->{sequence},
                $var2_fwd
            ) . "\n";
        }
    }
}

# Subroutines
sub read_file {
    my ($filename) = @_;
    my %data;
    
    open(my $fh, '<', $filename) or die "Cannot open $filename: $!\n";
    
    while (my $line = <$fh>) {
        chomp $line;
        next if $line =~ /^\s*$/;  # Skip empty lines
        
        my @fields = split /\s+/, $line;
        next if @fields < 3;  # Skip lines with insufficient columns
        
        my $marker_id = $fields[0];
        my $orientation = $fields[1];  # Second field contains orientation
        my $sequence = $fields[2];     # Third field contains sequence
        
        $data{$marker_id} = {
            orient   => $orientation,
            sequence => $sequence
        };
    }
    
    close($fh);
    return %data;
}

sub reverse_complement {
    my ($sequence) = @_;
    
    # Convert to uppercase for consistent processing
    $sequence = uc($sequence);
    
    # Reverse the sequence
    my $reversed = reverse($sequence);
    
    # Complement the bases
    $reversed =~ tr/ATCG/TAGC/;
    
    return $reversed;
}

sub print_usage {
    print <<EOF;
Usage: $0 <file1> <file2>

Compare two marker files and generate a comparison report.

Arguments:
    file1    First input file (marker_from file)
    file2    Second input file (marker_to file)

Input file format: markerID orientation sequence
- markerID: unique identifier (column 1)
- orientation: + or - (column 2)
- sequence: DNA sequence (column 3)

Output format: markerID compare len1 len2 orient var1 var2 var2_fwd
- orient field inherits orientation from file2
- var2_fwd: var2 sequence, reverse-complemented if orient is '-'

Example:
    $0 lis.marker_from lis.marker_to > report.tsv
EOF
}

