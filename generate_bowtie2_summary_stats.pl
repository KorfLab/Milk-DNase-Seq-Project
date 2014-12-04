#!/usr/bin/perl
#
# generate_bowtie2_summary_stats.pl
#
# A script to process output of bowtie2 on multiple FASTQ files to report on some
# summary stats of reads mapped to mm9 and mm10 versions of mouse genome.
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';

# one big hash to store everything (for mm9 and mm10)
my %bowtie2_stats;


# read all err files into list and loop over them

foreach my $file (glob("*.err")){
	
	warn "Processing $file\n";
	# extract prefix part of file name to get sample ID
	my ($id, $dnase) = $file =~ m/([A-Za-z0-9]+)_(\d+)/;
	my ($mouse_ver)  = $file  =~ m/(mm\d+)/;

	
	# now open err file and grab details to add to hash
	open(my $in, "<", $file) or die "Can't open $file";

	# only six lines per file, want info from four of them
	my $line1 = <$in>;
	my $line2 = <$in>;
	my $line3 = <$in>;
	my $line4 = <$in>;
	my $line5 = <$in>;
	
	my ($read_count)   = $line1 =~ m/(\d+) reads/;
	my ($unaligned)    = $line3 =~ m/(\d+) \(/;
	my ($aligned_once) = $line4 =~ m/(\d+) \(/;
	my ($aligned_many) = $line5 =~ m/(\d+) \(/;
		
	$bowtie2_stats{$id}{$dnase}{$mouse_ver}{read_count}   += $read_count;
	$bowtie2_stats{$id}{$dnase}{$mouse_ver}{unaligned}    += $unaligned;	
	$bowtie2_stats{$id}{$dnase}{$mouse_ver}{aligned_once} += $aligned_once;
	$bowtie2_stats{$id}{$dnase}{$mouse_ver}{aligned_many} += $aligned_many;
	close($in);
	
	# now process equivalent sam file to get to MAPQ scores
	my $sam_file = $file;
	$sam_file =~ s/err/sam/;
	warn "Processing $sam_file\n";
	open(my $sam, "<", $sam_file) or die "Can't open $sam_file";
	while(my $line = <$sam>){
		my @f = split(/\s+/, $line);
		my $mapq = $f[4];
		# add to main hash based on MAPQ score
		# categories taken from http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
		# may only want to keep reads with max score of 42
		if ($mapq == 42){
			$bowtie2_stats{$id}{$dnase}{$mouse_ver}{mapq_42}++;
		} elsif ($mapq >= 40){
			$bowtie2_stats{$id}{$dnase}{$mouse_ver}{mapq_40_41}++;		
		} elsif ($mapq >= 23){
			$bowtie2_stats{$id}{$dnase}{$mouse_ver}{mapq_23_39}++;		
		} elsif ($mapq >= 3){
			$bowtie2_stats{$id}{$dnase}{$mouse_ver}{mapq_3_22}++;		
		} elsif ($mapq >= 2){
			$bowtie2_stats{$id}{$dnase}{$mouse_ver}{mapq_2}++;		
		}else{
			$bowtie2_stats{$id}{$dnase}{$mouse_ver}{mapq_1}++;		
		}
	}
	close($sam);
}

# report on stats pooled for each sample for each genome version
print "ID\tDNAse_concentration\t";
print "Number_of_reads_mm9\t%_aligned_mm9\t";
print "%_MAPQ=42_mm9\t$%_MAPQ=1_mm9\t";

print "Number_of_reads_mm10\t%_aligned_mm10\t";
print "%_MAPQ=42_mm10\t$%_MAPQ=1_mm10\t";

print "\n";

foreach my $id (sort keys %bowtie2_stats){
	
	foreach my $dnase (sort {$a <=> $b} keys %{$bowtie2_stats{$id}}){
	
		print "$id\t$dnase\t";
		print "$bowtie2_stats{$id}{$dnase}{mm9}{read_count}\t";
		my $percent_aligned = sprintf("%.2f", ($bowtie2_stats{$id}{$dnase}{mm9}{aligned_once} + $bowtie2_stats{$id}{$dnase}{mm9}{aligned_many}) / $bowtie2_stats{$id}{$dnase}{mm9}{read_count} * 100);
		print "$percent_aligned\t";
	
		# print mapping percentages for some MAPQ score categories
		my $percent_42 = sprintf("%.2f", $bowtie2_stats{$id}{$dnase}{mm9}{mapq_42} / $bowtie2_stats{$id}{$dnase}{mm9}{read_count} * 100);
		print "$percent_42\t";

		my $percent_1 = sprintf("%.2f", $bowtie2_stats{$id}{$dnase}{mm9}{mapq_1} / $bowtie2_stats{$id}{$dnase}{mm9}{read_count} * 100);
		print "$percent_1\t";


		print "$bowtie2_stats{$id}{$dnase}{mm10}{read_count}\t";
		$percent_aligned = sprintf("%.2f", ($bowtie2_stats{$id}{$dnase}{mm10}{aligned_once} + $bowtie2_stats{$id}{$dnase}{mm10}{aligned_many}) / $bowtie2_stats{$id}{$dnase}{mm10}{read_count} * 100);
		print "$percent_aligned\t";


		# print mapping percentages for some MAPQ score categories
		$percent_42 = sprintf("%.2f", $bowtie2_stats{$id}{$dnase}{mm10}{mapq_42} / $bowtie2_stats{$id}{$dnase}{mm10}{read_count} * 100);
		print "$percent_42\t";

		$percent_1 = sprintf("%.2f", $bowtie2_stats{$id}{$dnase}{mm10}{mapq_1} / $bowtie2_stats{$id}{$dnase}{mm10}{read_count} * 100);
		print "$percent_1\t";

		print "\n";
	}
}