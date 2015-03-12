#!/usr/bin/perl
#
# generate_roc_data.pl
#
# A script to generate false positive and true positive data for an ROC curve
# True positive data based on overlap between DNase-seq data and STAT5 Chip-Seq peaks
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';

die "Usage:$0 <peak data file> <sam file>" unless (@ARGV == 2);

my ($peak_data_file, $sam_file) = @ARGV;

# first need to read peaks from STAT5
my %peaks;
my $peak_length;
read_stat5_peak_data($peak_data_file);

# want to assess things at various thresholds
foreach my $threshold (qw(5 10 15 20 25)){

	warn "Using threshold $threshold\n";
	# open SAM file
	open(my $in, "<", $sam_file) or die "Can't open $sam_file\n";
	while (my $line = <$in>){

		my @f = split(/\t/, $line);
		my ($chr, $beg, $seq) = ($f[2], $f[3], $f[9]);
		
		print "$chr\t$beg\t$seq\n";
		
	}
	close($in);
}


####################################
#
# S u b r o u t i n e s
#
####################################

sub read_stat5_peak_data{

	my ($file) = @_;

	die "$file should be a BED format file\n" if ($file !~ m/\.bed$/);

	# loop over BED file
	open(my $in, "<", $file) or die "Can't open $file\n";
	while (my $line = <$in>){

		next if ($line =~ m/^#/);
		next if ($line =~ m/^track/);
		my ($chr, $beg, $end, $score) = split(/\t/, $line);
		$peak_length = $end - $beg + 1;

		# add to hash of arrays
		push (@{${peaks}{$chr}}, $beg);
	}
	close($in);

	# now want to sort peaks on each chromosome for moe efficient look up later on
	foreach my $chr (sort keys %peaks){
		print "Sorting $chr peaks\n";
		my @sorted_coords = sort {$a <=> $b} @{$peaks{$chr}};
		@{$peaks{$chr}} = @sorted_coords;
	}

	foreach my $chr (sort keys %peaks){
		print "*$chr*\n@{$peaks{$chr}}\n\n";
	}
}