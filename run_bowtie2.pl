#!/usr/bin/perl
#
# run_bowtie2.pl
#
# A script to submit multiple bowtie jobs to cluster using qsub 
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';

die "Usage: $0 <number of cores>\n" unless (@ARGV == 1);

my ($threads) =  @ARGV;

my $mm9_index = "/share/tamu/Data/Bowtie2_indexes/mm9_index";
my $mm10_index = "/share/tamu/Data/Bowtie2_indexes/mm10_index";

foreach my $input_fastq_file (glob("*processed.fastq")){
	my $output_sam_file =  $input_fastq_file;
	my $error_file      =  $input_fastq_file;

	my $output_sam_file_mm9 =  $output_sam_file;
	my $output_sam_file_mm10 =  $output_sam_file;

	my $error_file_mm9  =  $error_file;
	my $error_file_mm10 =  $error_file;

	$output_sam_file_mm9    =~ s/_processed.fastq/_mm9.sam/;
	$output_sam_file_mm10   =~ s/_processed.fastq/_mm10.sam/;
	$error_file_mm9         =~ s/_processed.fastq/_mm9.err/;
	$error_file_mm10         =~ s/_processed.fastq/_mm10.err/;

	# for each FASTQ file we want to run against mm9 and mm10
	my $bowtie_command_mm9  = "bowtie2 -x $mm9_index -U $input_fastq_file -p $threads --very-sensitive --no-unal --no-hd -S $output_sam_file_mm9";
	my $bowtie_command_mm10 = "bowtie2 -x $mm10_index -U $input_fastq_file -p $threads --very-sensitive --no-unal --no-hd -S $output_sam_file_mm10";
	my $qsub_command_mm9  = "qsub -e $error_file_mm9 -V -l h_vmem=5G -pe make $threads -b y -cwd $bowtie_command_mm9 &";
	my $qsub_command_mm10 = "qsub -e $error_file_mm10 -V -l h_vmem=5G -pe make $threads -b y -cwd $bowtie_command_mm10 &";

	print "Running: $qsub_command_mm9\n\n";
	system($qsub_command_mm9) && die "Can't run $qsub_command_mm9\n";

	print "Running: $qsub_command_mm10\n\n";
	system($qsub_command_mm10) && die "Can't run $qsub_command_mm10\n";


}

exit;



