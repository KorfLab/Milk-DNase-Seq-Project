#!/usr/bin/perl
#
# sam_filter.pl
#
# Quick script to just extract entries from multiple sam files that match one chromosome
# and write to new output file
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';

die "Usage: $0 <chromosome to extract> <directory containing SAM files>\n" unless @ARGV == 2;

my ($chr, $sam_dir) = @ARGV;


# process the raw fastq files first
foreach my $file (glob("$sam_dir/*.sam")){
	warn "Processing $file\n";
	
	# form suitable output file name
	my $output = $file;
	$output =~ s/\.sam//;
	$output .= "_${chr}.sam";

	
	open(my $in,  "<", $file)   or die "Can't open $file\n";	
	open(my $out, ">", $output) or die "Can't write to $output\n";

	# process sam file, 
 	while (my $line = <$in>){
		my @f = split(/\s+/, $line);
		print $out "$line" if ($f[2] eq $chr);
    }
    close($in);
}


exit;

