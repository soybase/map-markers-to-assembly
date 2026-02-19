#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::File;

# Global variables
my ($vcf_file, $bed_file, $output_file, $reference_file, $help);

# Parse command line options
GetOptions(
    'vcf=s'       => \$vcf_file,
    'bed=s'       => \$bed_file,
    'output=s'    => \$output_file,
    'reference=s' => \$reference_file,
    'help|h'      => \$help
) or die "Error in command line arguments\n";

if ($help || !$vcf_file || !$bed_file || !$output_file || !$reference_file) {
    print_usage();
    exit;
}

# Main execution
main();

sub main {
    # Read genome assembly to get contig lengths
    print "Reading genome assembly file...\n";
    my $contig_lengths = read_genome_assembly($reference_file);
    
    # Read BED file data
    print "Reading BED file...\n";
    my $bed_data = read_bed_file($bed_file);
    
    # Process VCF file
    print "Processing VCF file...\n";
    process_vcf($vcf_file, $bed_data, $contig_lengths, $output_file);
    
    print "Processing complete. Output written to: $output_file\n";
}

sub read_genome_assembly {
    my ($reference_file) = @_;
    my %contig_lengths;

    # Open reference file (handle compression)
    my $ref_fh;
    if ($reference_file =~ /\.gz$/) {
        # Use MultiStream => 1 to handle bgzip files correctly
        $ref_fh = IO::Uncompress::Gunzip->new($reference_file, MultiStream => 1)
            or die "Cannot open compressed reference file '$reference_file': $GunzipError\n";
    } else {
        open $ref_fh, '<', $reference_file or die "Cannot open reference file '$reference_file': $!\n";
    }

    my $current_contig = '';
    my $current_length = 0;

    while (my $line = <$ref_fh>) {
        chomp $line;

        if ($line =~ /^>(.+)/) {
            # New contig header
            if ($current_contig && $current_length > 0) {
                $contig_lengths{$current_contig} = $current_length;
            }

            # Extract contig name (first word after >)
            my $header = $1;
            ($current_contig) = split /\s+/, $header;
            $current_length = 0;
        } elsif ($line !~ /^>/ && $line =~ /\S/) {
            # Sequence line - count nucleotides
            $line =~ s/\s//g;  # Remove any whitespace
            $current_length += length($line);
        }
    }

    # Don't forget the last contig
    if ($current_contig && $current_length > 0) {
        $contig_lengths{$current_contig} = $current_length;
    }

    close $ref_fh;

    print "Read " . scalar(keys %contig_lengths) . " contigs from genome assembly\n";
    return \%contig_lengths;
}

sub read_bed_file {
    my ($bed_file) = @_;
    my %bed_data;
    
    open my $bed_fh, '<', $bed_file or die "Cannot open BED file '$bed_file': $!\n";
    
    while (my $line = <$bed_fh>) {
        chomp $line;
        next if $line =~ /^#/ || $line =~ /^\s*$/;  # Skip comments and empty lines
        
        my @fields = split /\t/, $line;
        if (@fields < 7) {
            warn "Warning: BED line has fewer than 7 fields, skipping: $line\n";
            next;
        }
        
        my ($molecule, $start, $end, $marker_name, $score, $orientation, $allele) = @fields;
        
        # Store BED data keyed by marker name
        $bed_data{$marker_name} = {
            molecule    => $molecule,
            start       => $start,
            end         => $end,
            score       => $score,
            orientation => $orientation,
            allele      => $allele
        };
    }
    
    close $bed_fh;
    
    print "Read " . scalar(keys %bed_data) . " markers from BED file\n";
    return \%bed_data;
}

sub process_vcf {
    my ($vcf_file, $bed_data, $contig_lengths, $output_file) = @_;
    
    # Open VCF file (handle compression)
    my $vcf_fh;
    if ($vcf_file =~ /\.gz$/) {
        # Use MultiStream => 1 for bgzip compatibility
        $vcf_fh = IO::Uncompress::Gunzip->new($vcf_file, MultiStream => 1)
            or die "Cannot open compressed VCF file '$vcf_file': $GunzipError\n";
    } else {
        open $vcf_fh, '<', $vcf_file or die "Cannot open VCF file '$vcf_file': $!\n";
    }
    
    # First pass: collect all variant data for sorting
    my @header_lines;
    my @variant_records;
    my ($processed_variants, $matched_variants, $unmatched_variants) = (0, 0, 0);
    
    while (my $line = <$vcf_fh>) {
        chomp $line;
        
        # Collect header lines
        if ($line =~ /^#/) {
            push @header_lines, $line;
            next;
        }
        
        # Process variant lines
        my @vcf_fields = split /\t/, $line;
        if (@vcf_fields < 8) {
            warn "Warning: VCF line has fewer than 8 fields, skipping: $line\n";
            next;
        }
        
        $processed_variants++;
        
        my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, @format_fields) = @vcf_fields;
        
        # Use ID field as marker name to match with BED file
        my $marker_name = $id;
        
        if ($marker_name eq '.' || !exists $bed_data->{$marker_name}) {
            # No matching marker in BED file
            $unmatched_variants++;
            next;  # Skip this variant
        }
        
        $matched_variants++;
        
        # Get BED file information
        my $bed_info = $bed_data->{$marker_name};
        my $new_chrom = $bed_info->{molecule};
        my $new_pos = $bed_info->{start} + 1;  # Convert 0-based BED to 1-based VCF
        my $bed_allele = $bed_info->{allele};
        
        # Determine relationship between VCF ref and BED allele
        my ($new_ref, $new_alt) = process_alleles($ref, $alt, $bed_allele);
        
        # Store variant record for sorting
        push @variant_records, {
            chrom => $new_chrom,
            pos   => $new_pos,
            line  => join("\t", $new_chrom, $new_pos, $id, $new_ref, $new_alt, $qual, $filter, $info, @format_fields)
        };
    }
    
    close $vcf_fh;
    
    # Sort variants by chromosome and position
    @variant_records = sort {
        $a->{chrom} cmp $b->{chrom} || $a->{pos} <=> $b->{pos}
    } @variant_records;
    
    # Write output file
    open my $out_fh, '>', $output_file or die "Cannot create output file '$output_file': $!\n";
    
    # Write modified header
    write_modified_header($out_fh, \@header_lines, $contig_lengths);
    
    # Write sorted variant records
    for my $record (@variant_records) {
        print $out_fh $record->{line} . "\n";
    }
    
    close $out_fh;
    
    print "Summary:\n";
    print "  Total variants processed: $processed_variants\n";
    print "  Variants matched with BED: $matched_variants\n";
    print "  Variants without BED match: $unmatched_variants\n";
    print "  Variants written to output: " . scalar(@variant_records) . "\n";
}

sub write_modified_header {
    my ($out_fh, $header_lines, $contig_lengths) = @_;
    
    my $reference_written = 0;
    my $contig_written = 0;
    
    for my $line (@$header_lines) {
        # Replace reference line
        if ($line =~ /^##reference=/) {
            print $out_fh "##reference=$reference_file\n";
            $reference_written = 1;
        }
        # Skip existing contig lines - we'll write new ones
        elsif ($line =~ /^##contig=/) {
            # Skip existing contig lines
            next;
        }
        # Write contig lines after fileformat but before other headers
        elsif ($line =~ /^##fileformat=/ || ($line =~ /^##/ && !$contig_written && $reference_written)) {
            print $out_fh "$line\n";
            
            # Write contig lines after fileformat and reference
            if ($reference_written && !$contig_written) {
                write_contig_lines($out_fh, $contig_lengths);
                $contig_written = 1;
            }
        }
        else {
            print $out_fh "$line\n";
        }
    }
    
    # Add reference line if it wasn't in original header
    if (!$reference_written) {
        print $out_fh "##reference=$reference_file\n";
    }
    
    # Add contig lines if not written yet
    if (!$contig_written) {
        write_contig_lines($out_fh, $contig_lengths);
    }
}

sub write_contig_lines {
    my ($out_fh, $contig_lengths) = @_;
    
    my @sorted_contigs = sort keys %$contig_lengths;
    print "Writing " . scalar(@sorted_contigs) . " contig lines to header\n";
    
    # Sort contigs for consistent output
    for my $contig_id (@sorted_contigs) {
        my $length = $contig_lengths->{$contig_id};
        print $out_fh "##contig=<ID=$contig_id,length=$length>\n";
        print "  Written: ##contig=<ID=$contig_id,length=$length>\n" if @sorted_contigs <= 5; # Debug small cases
    }
}

sub process_alleles {
    my ($vcf_ref, $vcf_alt, $bed_allele) = @_;
    
    my $new_ref = $bed_allele;
    my $new_alt;
    
    # Check if VCF ref matches BED allele exactly
    if (uc($vcf_ref) eq uc($bed_allele)) {
        $new_alt = $vcf_alt;
    }
    # Check if VCF ref is reverse complement of BED allele
    elsif (uc($vcf_ref) eq uc(reverse_complement($bed_allele))) {
        $new_alt = reverse_complement_alt($vcf_alt);
    }
    # No simple relationship
    else {
        $new_alt = '.';
    }
    
    return ($new_ref, $new_alt);
}

sub reverse_complement {
    my ($seq) = @_;
    
    # Handle empty or undefined sequences
    return $seq unless defined $seq && $seq ne '';
    
    # Create reverse complement
    my $rev_comp = reverse($seq);
    $rev_comp =~ tr/ATCGatcg/TAGCtagc/;
    
    return $rev_comp;
}

sub reverse_complement_alt {
    my ($alt_field) = @_;
    
    # Handle missing alternate allele
    return '.' if !defined $alt_field || $alt_field eq '.';
    
    # Handle multiple alternate alleles separated by commas
    my @alts = split /,/, $alt_field;
    my @rev_comp_alts;
    
    for my $alt (@alts) {
        if ($alt eq '.') {
            push @rev_comp_alts, '.';
        } else {
            push @rev_comp_alts, reverse_complement($alt);
        }
    }
    
    return join(',', @rev_comp_alts);
}

sub print_usage {
    print <<EOF;
Usage: $0 --vcf <vcf_file> --bed <bed_file> --reference <reference_file> --output <output_file>

This script performs coordinate liftover between genome assemblies using:
  - A VCF file with variants in one assembly
  - A modified BED file with coordinates in another assembly
  - A reference genome assembly file for the target coordinates

Options:
  --vcf        Input VCF file (can be gzipped)
  --bed        Modified BED file (7 columns: molecule, start, end, marker-name, score, orientation, allele)
  --reference  Target genome assembly file (FASTA format, can be gzipped)
  --output     Output VCF file with updated coordinates
  --help       Show this help message

The script:
1. Updates the ##reference= header line with the provided reference file path
2. Calculates contig lengths from the reference and adds ##contig= header lines
3. Sorts output variants by chromosome/molecule and position
4. Matches variants by marker ID and updates coordinates and alleles accordingly

When the VCF reference allele is the reverse complement of the BED allele, 
the alternate allele is also reverse complemented.

Example:
  $0 --vcf variants.vcf.gz --bed markers.bed --reference target_genome.fasta.gz --output lifted_variants.vcf

EOF
}

