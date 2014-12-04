#!/usr/bin/perl
#
# filter_sam_files.pl
#
# A script to process bowtie 2 output of multiple SAM files and only keep 
# those with alignments above a certain mapping quality, and bin resulting SAM files
# based on sample ID
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';


die "Usage: $0 <MAPQ score threshold>\n" unless @ARGV == 1;

my ($mapq_threshold) = @ARGV;

# will need to create new SAM output files for each sample

my $file1 = "1b_5_mm10_MAPQ${mapq_threshold}.sam";
my $file2 = "1b_10_mm10_MAPQ${mapq_threshold}.sam";
my $file3 = "1b_50_mm10_MAPQ${mapq_threshold}.sam";
my $file4 = "2_10_mm10_MAPQ${mapq_threshold}.sam";
my $file5 = "2_20_mm10_MAPQ${mapq_threshold}.sam";
my $file6 = "2_50_mm10_MAPQ${mapq_threshold}.sam";
my $file7 = "8Lvr_10_mm10_MAPQ${mapq_threshold}.sam";
my $file8 = "8Lvr_20_mm10_MAPQ${mapq_threshold}.sam";
my $file9 = "8Lvr_50_mm10_MAPQ${mapq_threshold}.sam";
my $file10 = "8Mg_5_mm10_MAPQ${mapq_threshold}.sam";
my $file11 = "8Mg_10_mm10_MAPQ${mapq_threshold}.sam";
my $file12 = "8Mg_20_mm10_MAPQ${mapq_threshold}.sam";
my $file13 = "8Mg_50_mm10_MAPQ${mapq_threshold}.sam";
my $file14 = "16P_5_mm10_MAPQ${mapq_threshold}.sam";
my $file15 = "16P_10_mm10_MAPQ${mapq_threshold}.sam";
my $file16 = "16P_20_mm10_MAPQ${mapq_threshold}.sam";
my $file17 = "16P_50_mm10_MAPQ${mapq_threshold}.sam";

open(my $out1, ">", $file1) or die "Can't write to $file1";
open(my $out2, ">", $file2) or die "Can't write to $file2";
open(my $out3, ">", $file3) or die "Can't write to $file3";
open(my $out4, ">", $file4) or die "Can't write to $file4";
open(my $out5, ">", $file5) or die "Can't write to $file5";
open(my $out6, ">", $file6) or die "Can't write to $file6";
open(my $out7, ">", $file7) or die "Can't write to $file7";
open(my $out8, ">", $file8) or die "Can't write to $file8";
open(my $out9, ">", $file9) or die "Can't write to $file9";
open(my $out10, ">", $file10) or die "Can't write to $file10";
open(my $out11, ">", $file11) or die "Can't write to $file11";
open(my $out12, ">", $file12) or die "Can't write to $file12";
open(my $out13, ">", $file13) or die "Can't write to $file13";
open(my $out14, ">", $file14) or die "Can't write to $file14";
open(my $out15, ">", $file15) or die "Can't write to $file15";
open(my $out16, ">", $file16) or die "Can't write to $file16";
open(my $out17, ">", $file17) or die "Can't write to $file17";


# read all mm10 SAM files into list and loop over them

my $discarded_alignments = 0;
my $alignment_counter = 0;

foreach my $file (glob("*mm10.sam")){
	
	warn "Processing $file\n";
	# extract prefix part of file name to get sample ID
	my ($id, $dnase) = $file =~ m/([A-Za-z0-9]+)_(\d+)/;
	

	open(my $sam, "<", $file) or die "Can't open $file";
	SAM: while(my $line = <$sam>){
		$alignment_counter++;
		my @f = split(/\s+/, $line);
		my $mapq = $f[4];

		# skip if MAPQ score is too low
		if ($mapq < $mapq_threshold){
			$discarded_alignments++;
			next SAM;
		} 
		
		# now print to correct file handle
		
		if ($id eq '1b' and $dnase == 5){
			print $out1 "$line";
		} elsif($id eq '1b' and $dnase == 10){
			print $out2 "$line";	
		} elsif($id eq '1b' and $dnase == 50){
			print $out3 "$line";	
		} elsif($id eq '2' and $dnase == 10){
			print $out4 "$line";	
		} elsif($id eq '2' and $dnase == 20){
			print $out5 "$line";	
		} elsif($id eq '2' and $dnase == 50){
			print $out6 "$line";	
		} elsif($id eq '8Lvr' and $dnase == 10){
			print $out7 "$line";	
		} elsif($id eq '8Lvr' and $dnase == 20){
			print $out8 "$line";	
		} elsif($id eq '8Lvr' and $dnase == 50){
			print $out9 "$line";	
		} elsif($id eq '8Mg' and $dnase == 5){
			print $out10 "$line";	
		} elsif($id eq '8Mg' and $dnase == 10){
			print $out11 "$line";	
		} elsif($id eq '8Mg' and $dnase == 20){
			print $out12 "$line";	
		} elsif($id eq '8Mg' and $dnase == 50){
			print $out13 "$line";	
		} elsif($id eq '16P' and $dnase == 5){
			print $out14 "$line";	
		} elsif($id eq '16P' and $dnase == 10){
			print $out15 "$line";	
		} elsif($id eq '16P' and $dnase == 20){
			print $out16 "$line";	
		} elsif($id eq '16P' and $dnase == 50){
			print $out17 "$line";	
		} else {
			die "Shouldn't get here!\n$line\n";
		}
		
	}
	close($sam);
}

warn "Discarded $discarded_alignments out of $alignment_counter alignments with MAPQ scores < $mapq_threshold\n";



close($out1);
close($out2);
close($out3);
close($out4);
close($out5);
close($out6);
close($out7);
close($out8);
close($out9);
close($out10);
close($out11);
close($out12);
close($out13);
close($out14);
close($out15);
close($out16);
close($out17);

exit;