#!/usr/bin/perl
#
# chromosome_sizes.pl
#
# want to find out how big each gzipped chromosome file is (plus number of N bases)
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';
use FAlite;

die "Usage: $0 <directory of gzipped FASTA files>\n" unless @ARGV == 1;

my ($fasta_dir) = @ARGV;

# mm9:  Total Bases in Assembly 2,745,142,291 Total Non-N Bases in Assembly 2,648,522,751
# From my script: 2,725,765,481 bp

# mm10: Total Bases in Assembly 2,798,785,524 Total Non-N Bases in Assembly 2,719,482,043
# From my script 2,730,871,774


foreach my $file (glob("$ARGV[0]/*.fa.gz")){
	warn "Processing $file\n";
	my $chromosome_size = 0;
	open(my $FASTA, "gunzip -c $file | ") or die "Can't open pipe: $? $!\n";
	my $FA = new FAlite ($FASTA);
	my $n_count;
    while (my $entry = $FA->nextEntry) {
    	$chromosome_size += length($entry->seq);
    	my $n = ($entry->seq =~ tr/N/N/);
    	$n_count += $n;
    }
   	close($FASTA);
	print "$file\t$chromosome_size\t$n_count\n";
}

exit;