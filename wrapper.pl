#!/usr/bin/perl
#
# wrapper.pl
#
# Simple script to run read_depth.pl script across a bunch of samples
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';

die "Usage: $0 <bin size>\n" unless @ARGV == 1;

my ($bin_size) = @ARGV;

my $script = "read_depth.pl";

foreach my $file (glob("*MAPQ24.sam")){
        
    warn "\nProcessing $file\n";

	# extract prefix part of file name and the sample ID
	my ($file_prefix) = $file =~ m/(.*MAPQ24)/;
	
	my $output_file = $file_prefix . "_bin${bin_size}.tsv";
	my $command = "$script $bin_size 0 $file > $output_file";
	warn "$command\n";
	system ($command) && die "Can't run $command\n";
}

warn "Finished\n";
exit;