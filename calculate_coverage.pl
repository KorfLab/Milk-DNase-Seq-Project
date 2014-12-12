#!/usr/bin/perl
#
# calculate_coverage.pl
#
# want to find out how size of alignments (from SAM files) and size of pre- and post-
# processed FASTQ files, compares to sizes of genome (mm9 and mm10)
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';
use FAlite;

die "Usage: $0 <directory of SAM files> <directory of FASTQ files>\n" unless @ARGV == 3;

my ($sam_dir, $fastq_dir, $seqs) = @ARGV;

my $genome_size_mm9  = 2725765481;
my $genome_size_mm10 = 2730871774;
my $mapq_threshold = 24;

my %sequence_depth;


# process the raw fastq files first
FILE: foreach my $file (glob("$fastq_dir/*.fastq.gz")){
	#next unless ($file =~ m/8Mg_50/);
	warn "Processing $file\n";
	
	my $fastq_seq_counter = 0;

	# extract prefix part of file name and the sample ID
	my ($file_prefix) = $file =~ m/([A-Za-z0-9]+_\d+_[ACGT]*_L\d+_R1_\d+)/;
	my ($id)          = $file =~ m/([A-Za-z0-9]+_\d+)/;

	# is this a raw FASTQ file, or one processed (by Scythe and Sickle)?
	my $type = 'raw';
	$type = 'processed' if ($file =~ m/processed/);

	# process fastq file, just want to know length of line 2 (sequence)
	open(my $FASTQ, "gunzip -c $file | ") or die "Can't open pipe: $? $!\n";
 	FASTQ: while (my $line1 = <$FASTQ>){
 		my $line2 = <$FASTQ>;
		my $line3 = <$FASTQ>;
		my $line4 = <$FASTQ>;
		chomp($line2);

		$fastq_seq_counter++;
		last FASTQ if ($fastq_seq_counter > $seqs);

		# add length to hash		
		$sequence_depth{$id}{$type} += length($line2);			
		
		
    }
    close($FASTQ);

    # now open equivalent SAM files (but only do this for each raw FASTQ file)
    # don't want to process files twice 
    next FILE if ($type eq 'processed');
    my $sam_file_mm9  = $file_prefix . "_mm9.sam.gz";
	my $sam_file_mm10 = $file_prefix . "_mm10.sam.gz";

	foreach my $sam_file ($sam_file_mm9, $sam_file_mm10){
		warn "Processing $sam_file\n";
		process_sam_file($sam_file, $id);
	}

}

# level of precision for coverage values
my $precision = 5;

# print out final coverage values for raw and processed sequences
# and for mm9 and mm10
# and for various levels of SAM filtering (only high MAPQ scores and non-ChrM mappings)

print "ID\t";
print "Raw_coverage_mm9\tProcessed_coverage_mm9\tSAM_mm9\tSAM_mm9_MAPQ\tSAM_mm9_non_M\t";
print "Raw_coverage_mm10\tProcessed_coverage_mm10\tSAM_mm10\tSAM_mm10_MAPQ\tSAM_mm10_non_M\n";

foreach my $id (sort keys %sequence_depth){
	my $raw_depth            = $sequence_depth{$id}{raw};
	my $processed_depth      = $sequence_depth{$id}{processed};
	my $raw_sam_depth_mm9    = $sequence_depth{$id}{sam_raw}{mm9};
	my $raw_sam_depth_mm10   = $sequence_depth{$id}{sam_raw}{mm10};
	my $sam_depth_mapq_mm9   = $sequence_depth{$id}{sam_high_mapq}{mm9};
	my $sam_depth_mapq_mm10  = $sequence_depth{$id}{sam_high_mapq}{mm10};
	my $sam_depth_non_M_mm9  = $sequence_depth{$id}{sam_high_mapq_non_chrM}{mm9};
	my $sam_depth_non_M_mm10 = $sequence_depth{$id}{sam_high_mapq_non_chrM}{mm10};

	# print "$raw_depth\t$processed_depth\t$raw_sam_depth_mm9\t$sam_depth_mapq_mm9\t$sam_depth_non_M_mm9\n";
	# print "$raw_depth\t$processed_depth\t$raw_sam_depth_mm10\t$sam_depth_mapq_mm10\t$sam_depth_non_M_mm10\n";
	# exit;

	my $raw_coverage_mm9         = sprintf("%.${precision}f", $raw_depth / $genome_size_mm9); 
	my $processed_coverage_mm9   = sprintf("%.${precision}f", $processed_depth / $genome_size_mm9); 
	my $raw_sam_coverage_mm9     = sprintf("%.${precision}f", $raw_sam_depth_mm9 / $genome_size_mm9); 
	my $sam_mapq_coverage_mm9    = sprintf("%.${precision}f", $sam_depth_mapq_mm9 / $genome_size_mm9); 
	my $sam_non_M_coverage_mm9   = sprintf("%.${precision}f", $sam_depth_non_M_mm9 / $genome_size_mm9); 

	my $raw_coverage_mm10        = sprintf("%.${precision}f", $raw_depth / $genome_size_mm10); 
	my $processed_coverage_mm10  = sprintf("%.${precision}f", $processed_depth / $genome_size_mm10); 
	my $raw_sam_coverage_mm10    = sprintf("%.${precision}f", $raw_sam_depth_mm10 / $genome_size_mm10); 
	my $sam_mapq_coverage_mm10   = sprintf("%.${precision}f", $sam_depth_mapq_mm10 / $genome_size_mm10); 
	my $sam_non_M_coverage_mm10  = sprintf("%.${precision}f", $sam_depth_non_M_mm10 / $genome_size_mm10); 

	# final output
	print "$id\t";
	print "$raw_coverage_mm9\t$processed_coverage_mm9\t";
	print "$raw_sam_coverage_mm9\t$sam_mapq_coverage_mm9\t$sam_non_M_coverage_mm9\t";
	print "$raw_coverage_mm10\t$processed_coverage_mm10\t";
	print "$raw_sam_coverage_mm10\t$sam_mapq_coverage_mm10\t$sam_non_M_coverage_mm10\t";
	print "\n";
}

exit;


### subroutines
sub process_sam_file{
	my ($sam_file, $id) = @_;
	
	my $line_counter = 0;

	my $mouse_ver;
	$mouse_ver = 'mm9'  if ($sam_file =~ m/mm9/);
	$mouse_ver = 'mm10' if ($sam_file =~ m/mm10/);
	open(my $sam, "gunzip -c $sam_file | ") or die "Can't open pipe: $? $!\n";
	
	SAM: while(my $line = <$sam>){
		$line_counter++;
		last SAM if ($line_counter > $seqs);
			
	 	my @f = split(/\s+/, $line);
		my ($chr, $mapq, $sequence) = ($f[2], $f[4], $f[9]);
		my $seq_length = length($sequence);
		
		$sequence_depth{$id}{sam_raw}{$mouse_ver} += $seq_length;
		
		# also store results for high quality mappings
		# and then non-ChrM matches
		if ($mapq > $mapq_threshold){
			$sequence_depth{$id}{sam_high_mapq}{$mouse_ver} += $seq_length;

			if ($chr ne 'chrM'){
				$sequence_depth{$id}{sam_high_mapq_non_chrM}{$mouse_ver} += $seq_length;
			}
		}		
	}
	close($sam);
}
