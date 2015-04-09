#!/usr/bin/perl
#
# wrapper3.pl
#
# Simple script to run TopHat shell script across a bunch of samples
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';

die "Usage: $0 <pattern to match file names>" unless (@ARGV == 1);

my ($pattern) = @ARGV;

foreach my $file_A (glob("$pattern*1_sequence_processed.fastq")){
        
    warn "\nProcessing $file_A\n";

    # make file name for matching pair
	my $file_B = $file_A;
	$file_B =~ s/1_sequence/2_sequence/;

	# set up TopHat output directory name
	my $tophat_output = "/share/tamu/Analysis/TopHat_output/" . $file_A;
	$tophat_output =~ s/_1_sequence_processed.fastq//;

	my $command = "qsub -S /bin/bash -pe threaded 1 -l h_vmem=8G -M keith\@bradnam.co -m be -N krb_tophat_run /share/tamu/Code/run_tophat.sh ";
	$command .= " $file_A $file_B $tophat_output";
	warn "About to run $command\n";
	system ($command) && die "Can't run $command\n";
}

warn "Finished\n";
exit;