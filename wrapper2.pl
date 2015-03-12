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


foreach my $file_A (glob("*1_sequence.fastq.gz")){
        
    warn "\nProcessing $file_A\n";

    # make file name for matching pair
	my $file_B = $file_A;
	$file_B =~ s/1_sequence/2_sequence/;

	# extract prefix part of file name and the sample ID
	my ($file_prefix_A) = $file_A =~ m/(.*).fastq.gz/;
	my ($file_prefix_B) = $file_B =~ m/(.*).fastq.gz/;

	# set up scythe output file names
	my $scythe_output_A = $file_prefix_A . "_scythe.fastq";
	my $scythe_output_B = $file_prefix_B . "_scythe.fastq";

	# set up sickle output file names
	my $sickle_output_A = $file_prefix_A . "_processed.fastq";
	my $sickle_output_B = $file_prefix_B . "_processed.fastq";
	# deal with awkward singles file
	my $sickle_output_C = $file_prefix_A;
	$sickle_output_C =~ s/1_sequence/sequence_singles.fastq/;


	my $command = "qsub -S /bin/bash -pe threaded 2 -M keith\@bradnam.co -N krb_scythe_and_sickle_run /share/tamu/Code/run_scythe_and_sickle2.sh ";
	$command .= "$file_A $file_B $scythe_output_A $scythe_output_B $sickle_output_A $sickle_output_B $sickle_output_C";
	warn "$command\n";
	system ($command) && die "Can't run $command\n";
}

warn "Finished\n";
exit;