#!/usr/bin/perl
#
# join.pl
#
# A script to join a bunch of files together using Unix join command
#
# Author: Keith Bradnam, Genome Center, UC Davis
# This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
# This software is provided AS IS, without warranty of any kind.

use strict;
use warnings FATAL => 'all';

die "Usage: $0 <pattern for files to join>" unless @ARGV == 1;
my ($pattern) = @ARGV;

# get list of files to join
my @files = glob("$pattern");

# to start we need to join the first two files separately
my $file_A = shift(@files);
my $file_B = shift(@files);

# will also need to add header info into file (columns 'GeneID' and <filename>) if not present
# this will be done in separate files so as to keep the originals intact
my $file_A_tmp = "${file_A}.tmp";
my $file_B_tmp = "${file_B}.tmp";

add_header($file_A, $file_A_tmp);
add_header($file_B, $file_B_tmp);

# now join two files
my $tmp_join_file = "tmp_join_file.txt";

my $command = "join -t \$'\\t' $file_A_tmp $file_B_tmp > $tmp_join_file";
warn "Running $command\n";
system($command) && die "Can't run $command";

# clean up
unlink($file_A_tmp) or warn "ERROR: Can't remove $file_A_tmp\n";
unlink($file_B_tmp) or warn "ERROR: Can't remove $file_B_tmp\n";


# can now loop over remaining files and perform extra join commands
foreach my $file (@files){
	my $file_tmp = "${file}.tmp";
	add_header($file, $file_tmp);
	warn "Joining $file_tmp to $tmp_join_file\n";
	my $tmp_join_file2 = "tmp_join_file2.txt";
	my $command = "join -t \$'\\t' $tmp_join_file $file_tmp > $tmp_join_file2";
	warn "Running $command\n";
	system($command) && die "Can't run $command";

	# now overwrite first temp file
	$command = "mv $tmp_join_file2 $tmp_join_file";

	warn "Running $command\n";
	system($command) && die "Can't run $command";

	# clean up
	unlink($file_tmp) or warn "ERROR: Can't remove $file_tmp\n";


}

exit;

sub add_header{
	my ($infile, $outfile) = @_;

	# modify filename to remove .txt
	my $infile_truncated = $infile;
	$infile_truncated =~ s/\.txt//;
	
	# add header to new temp file
	my $command = "echo -e 	\"GeneID\t$infile_truncated\" > $outfile";
	warn "Running $command\n";
	system($command) && die "Can't run $command";

	# concatenate all of original file into new file
	$command = "cat $infile >> $outfile";
	warn "Running $command\n";
	system($command) && die "Can't run $command";
}

