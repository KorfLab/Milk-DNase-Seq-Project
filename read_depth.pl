#!/usr/bin/perl
#
# read_depth.pl
#
# Simple script to bin read counts from a SAM file
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';

die "Usage: $0 <bin size> <threshold> <file>\n" unless @ARGV == 3;

my ($bin_size, $threshold, $file) = @ARGV;

# need a hash of chromosome sizes?

# hash of arrays to store information
# primary key will be chromosome name
my %mapped_read_counts;	


open(my $in,  "<", $file)   or die "Can't open $file\n";	
	
# process sam file, 
while (my $line = <$in>){
	my @f = split(/\s+/, $line);
	my ($strand, $chr, $pos) = ($f[1], $f[2], $f[3]);	

	# bin index is provided by $pos/$bin_size
	my $index = $pos/$bin_size;
	my $floor = int($index);
	
	$mapped_read_counts{$chr}[$index]++; 
#	print "$chr\t$pos\t$index\t$floor\t$mapped_read_counts{$chr}[$index]\n";
}
close($in);


# want to track some stats. E.g. 
# 	% of bins that have zero counts
#   max read count for any bin in a sample
#   average read count per bin
my $all_bin_count = 0;
my $all_read_count = 0;
my $zero_bin_count = 0;
my $max_read_count = 0;
my $max_read_bin;

foreach my $chr (sort keys %mapped_read_counts){

	for (my $i = 0; $i < @{$mapped_read_counts{$chr}}; $i++){

		$all_bin_count++;

		if (not defined $mapped_read_counts{$chr}[$i]){
			$mapped_read_counts{$chr}[$i] = 0;
		}
		
		# extract basic stats of this bin
		my $beg = $i * $bin_size;
		my $end = $beg + $bin_size -1;
		my $count = $mapped_read_counts{$chr}[$i];
		$all_read_count += $count;
		$zero_bin_count++ if ($count == 0);

		if ($count > $max_read_count){
			$max_read_count = $count;
			$max_read_bin = "$chr $beg-$end";
		}
		
		print "$chr\t$beg\t$end\t$count\n" if ($count >= $threshold);
	}
	
}

my $average_bin_count = sprintf("%.2f", $all_read_count / $all_bin_count); 
my $zero_bin_percent  = sprintf("%.2f", $zero_bin_count / $all_bin_count * 100);

warn "Average read count per bin was $average_bin_count\n";
warn "Max read count observed was in bin $max_read_bin: $max_read_count reads\n";
warn "$zero_bin_percent%  of bins had zero read counts\n\n";
exit;

