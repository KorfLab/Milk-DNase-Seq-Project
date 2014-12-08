Milk-DNase-Seq-Project
======================

Analysis of milk DNase-Seq data in mouse and cow.


# Setup of directories for raw data #

In the shared folder, we will create:

	Data/
	Data/Genomes/
	Data/Genomes/Cow/
	Data/Genomes/Mouse/	
	Data/DNAse-Seq/
	Data/DNAse-Seq/Cow/
	Data/DNAse-Seq/Mouse/
	Data/RNA-Seq/
	Data/RNA-Seq/Cow/
	Data/RNA-Seq/Mouse/

	
`Genomes` subdirectory will ultimately contain reference genome sequences and annotations, whereas `DNAse-Seq` and `RNA-Seq` directories will contain our raw experimental data (in date-versioned subdirectories).
 
I've also added a `Data/RNA-Seq/Cow/Metadata` for a couple of files Danielle provided about the Cow RNA-Seq data:

`CowLactationRNASeqMetadata.txt` - A guess of what the sample IDs really mean
`CowLactationRNASeqMappingStats.tsv` - tab delimited file containing mapping statistics from Baylor

## DNAse-Seq data for mouse ##
The data exists as a set of zip files (one per sample). Each zip file contains several FASTQ files. Need to unzip everything (until we have enough secondary output files):

```bash
cd /share/tamu/Data/DNase-Seq/Mouse
# need to quote names for expansion to work with unzip (don't know why)
$ unzip '*.zip'
```


## Genome data ##
Genome datasets for mouse and cow were downloaded from the UCSC Genome Browser FTP site using rsync:

```bash
cd Data/Genomes/Mouse/mm9
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes .
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/mm9/vsBosTau7 .
cd ../mm10
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes .
cd ../../Cow/
mkdir bosTau7
cd bosTau7
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/bosTau7/bigZips ./

```

Note that two versions of the mouse genome were downloaded. Also, for mm9 I downloaded the comparison alignments to bosTau7. Initially, just the 'chromosomes' subdirectories of the mouse URLs were copied.

## Extracting adaptor sequences ##

From supplied Illumina PDF file, I made a single FASTA file containing sequence of the Universal TruSeq adapter, and the 24 TruSeq index adapters, and placed in  `Data/DNase-Seq/Mouse/adapter_sequences.fasta`

I then extracted the barcodes from the index adapters as follows:

```bash
cat adapter_sequences.fasta | sed 's/GATCGGAAGAGCACACGTCTGAACTCCAGTCAC//' | grep -v ">" | tail -n 24 | cut -c 1-6 > barcodes.txt
```


# Other directories for project #

A top level `Code` folder will contain scripts used to analyze data. This may be expanded to contain subdirectories (if we have some code in Git repositories).

A top level `Packages` folder will be used to install pre-existing bioinformatics tools. Until `git` is available on the cluster, may have to `git clone` locally and then `scp` files to our machine on the cluster.

A top level `bin` directory will contain symbolic links to executables in either `Code` or `Packages`. This `bin` directory should be added to users' $PATH.

A top level `Analysis` directory which will contain results from any steps we run.


# Making test dataset #

While testing Scythe and other tools, it will be good to just one work with a few FASTQ files. Make copy of `Sample_2_10.zip` sequences as this is the smallest file:

```bash
cd Analyis
mkdir Test
cp ../../Data/DNase-Seq/Mouse/Sample_2_10.zip .
unzip Sample_2_10.zip
cd Sample_2_10
```


## Checking barcodes in test dataset ##

Presumably barcodes in file names should match barcodes in the FASTQ identifiers in the respective files. But this is not what happens:

```bash
cd Analysis/Test
mkdir Barcode_check
cd Barcode_check
gunzip -c 2_10_AGTTCC*.fastq.gz | grep "@HWI" | sed 's/.* 1:N:0://' | sort | uniq -c > barcodes_in_identifiers.txt

cat barcodes_in_identifiers.txt

  5333 AATTCC
   5726 ACTTCC
   3433 AGATCC
   3424 AGCTCC
   6579 AGGTCC
    614 AGNTCC
   4143 AGTACC
   5102 AGTCCC
    713 AGTGCC
    172 AGTNCC
  34555 AGTTAC
  60077 AGTTCA
14783490 AGTTCC
  12199 AGTTCG
  21071 AGTTCT
   3341 AGTTGC
  23722 AGTTTC
  20326 ATTTCC
  20336 CGTTCC
   4689 GGTTCC
   2378 TGTTCC
   
```bash

Two of these contain an N character in the barcode. Are the non AGTTCC barcodes all part of the official set of barcodes that I extracted above?

```bash
cat barcodes_in_identifiers.txt | sed 's/ .* //' | grep -v AGTTCC > extra_barcodes.txt
grep -f extra_barcodes.txt /share/tamu/Data/DNase-Seq/Mouse/adapter_sequences.fasta


```

No matches at all! Though the two N cases would match if you treat them as ambiguous bases.

Main run:
```bash
scythe -a ../../../Data/DNase-Seq/Mouse/adapter_sequences.fasta -o 2_10_AGTTCC_L008_R1_001_trimmed.fastq 2_10_AGTTCC_L008_R1_001.fastq.gz

prior: 0.300

Adapter Trimming Complete
contaminated: 264003, uncontaminated: 3735997, total: 4000000
contamination rate: 0.066001

```

# Installing Scythe (v0.991) #
1. Cloned locally from https://github.com/vsbuffalo/scythe
2. `mkdir -p Packages/Scythe`
3. Use `scp` to copy local `scythe` directory to `Packages/Scythe/c128b19c65` # use SHA-1 hash for directory name
4. `cd Packages/Scythe`
5. `ln -s c128b19c65 current`
6. `cd current`
7. `make all`
8. No errors reported at this point, so…
9. `cd /share/tamu/bin`
10. `ln -s ../Packages/Scythe/current/scythe`


## Testing Scythe: part 1 ##
1. Using supplied sample files that are part of Scythe installation
2. `mkdir -p Analysis/Test/Scythe_test`
3. `cd Analysis/Test/Scythe_test`
4. Run Scythe:

```bash
scythe -a ../../../Packages/Scythe/current/illumina_adapters.fa -o trimmed_seqs.fastq ../../../Packages/Scythe/current/testing/reads.fastq
prior: 0.300

Adapter Trimming Complete
contaminated: 1329, uncontaminated: 8671, total: 10000
contamination rate: 0.132900
```

Resulting output file does indeed to be trimmed:

```bash
ls -lh trimmed_seqs.fasta ../../../Packages/Scythe/current/testing/reads.fastq
-rw-r--r--+ 1 keith dglemay 2.3M Oct 17 06:45 ../../../Packages/Scythe/current/testing/reads.fastq
-rw-r--r--+ 1 keith dglemay 2.0M Oct 17 07:02 trimmed_seqs.fastq
```


## Testing Scythe: part 2 ##

Now to try with some real DNAse-Seq data:

```bash
cd Analysis/Test/Scythe_test
ln -s ../Sample_2_10/2_10_AGTTCC_L008_R1_001.fastq.gz
ln -s ../../../Data/DNase-Seq/Mouse/adapter_sequences.fasta
scythe -a adapter_sequences.fasta -o 2_10_AGTTCC_L008_R1_001_trimmed.fastq 2_10_AGTTCC_L008_R1_001.fastq.gz

prior: 0.300

Adapter Trimming Complete
contaminated: 264003, uncontaminated: 3735997, total: 4000000
contamination rate: 0.066001


head -n 16 2_10_AGTTCC_L008_R1_001_trimmed.fastq
@HWI-ST976:178:c53v2acxx:8:1101:1273:2083 1:N:0:AGTTCC
AGGTTTCAGTTGATTGTCCTGTTTTCACTCACTGATTAGACAGTTAAATAA
+
@CCDFFFFHFHDBIIJGGGEGEHHJIJJIJJJJJJJJIBFHHIBGIIIJFI
@HWI-ST976:178:c53v2acxx:8:1101:1277:2127 1:N:0:AGTTCC
N
+
B
@HWI-ST976:178:c53v2acxx:8:1101:1317:2174 1:N:0:AGTTCC
TGACGACTTGAAAAATGACGAAATCACTAAAAAACGTGAAAAATGAGAAAT
+
@@@FDDDDFHFHHJGBGHIGHGGIJJJJJJJJJJIGBGIIIIIGGHGIJEG
@HWI-ST976:178:c53v2acxx:8:1101:1453:2223 1:N:0:AGTTCC
N
+
B
```

This suggests that when trimming occurs, the sequence is being replaced with a single N (subsequently confirmed with Vince Buffalo). 

Now to see how much sequence is being trimmed (using some awk skillz to just grab every second line from file):

```bash
awk 'NR %4 == 2 {print length($0); } ' 2_10_AGTTCC_L008_R1_001_trimmed.fastq | sort | uniq -c
 238554 1
   1124 36
   1658 37
   1618 38
   1546 39
   2332 40
   2087 41
   2179 42
   3796 43
   4059 44
   5050 45
3735997 51

```

So about 93% (3735997/(16000000/4) of sequences undergo no trimming, about 6% of sequences are completely removed by trimming, and about 0.5% have an intermediate number of bases removed.


# Installing Sickle (v1.33) #
1. Cloned locally from https://github.com/najoshi/sickle
2. `mkdir -p Packages/Sickle`
3. Use `scp` to copy local `sickle` directory to `Packages/Sickle/7667f14` # use SHA-1 hash for directory name
4. `cd Packages/Sickle`
5. `ln -s 7667f14 current`
6. `cd current`
7. `make`
8. No errors reported at this point, so…
9. `cd /share/tamu/bin`
10. `ln -s ../Packages/Sickle/current/sickle`

## Testing Sickle ##

This is single-end data, so we can use the 'se' mode of Sickle and test on a raw FASTQ file. Will use 'sanger' mode for quality encoding (which is equivalent to CASAVA >= 1.8).

```bash
cd /share/tamu/Analysis/Test
mkdir Sickle_test
cd Sickle_test
ln -s ../Sample_2_10/2_10_AGTTCC_L008_R1_001.fastq.gz

sickle se -f 2_10_AGTTCC_L008_R1_001.fastq.gz -o 2_10_AGTTCC_L008_R1_001_trimmed.fastq -t sanger

SE input file: 2_10_AGTTCC_L008_R1_001.fastq.gz

Total FastQ records: 4000000
FastQ records kept: 3968789
FastQ records discarded: 31211
```

Sickle removes sequences if the trimmed sequence is less than a specified length (default = 20 bp), so just under 1% of the sequences were completely removed for this file.


Now can try the same thing but on the already adapter-trimmed sequence from Scythe:

```bash
sickle se -f ../Scythe_test/2_10_AGTTCC_L008_R1_001_trimmed.fastq -o 2_10_AGTTCC_L008_R1_001_trimmedx2.fastq -t sanger

SE input file: ../Scythe_test/2_10_AGTTCC_L008_R1_001_trimmed.fastq

Total FastQ records: 4000000
FastQ records kept: 3739975
FastQ records discarded: 260025

```

Now ends up discarding about 7% of input reads. Finally, try the `-n` option to also trim at first appearance of an N:

```bash
sickle se -f ../Scythe_test/2_10_AGTTCC_L008_R1_001_trimmed.fastq -o 2_10_AGTTCC_L008_R1_001_trimmedx2_no_n.fastq -t sanger -n

SE input file: ../Scythe_test/2_10_AGTTCC_L008_R1_001_trimmed.fastq

Total FastQ records: 4000000
FastQ records kept: 3738889
FastQ records discarded: 261111
```

This discards 1,086 more sequences than before, suggesting that there must have been some sequences which have Ns, but which otherwise were not low quality (and are not the single-N sequences generated by Scythe):

```bash
awk 'NR % 4 == 2' ../Scythe_test/2_10_AGTTCC_L008_R1_001_trimmed.fastq  | grep N | grep -vE "^N$" | wc -l
1513
```

So, *after* trimming with Scythe, there were 1,513 sequences with at least one ambiguous base (exlcuding those with *only* one ambiguous base), but the earlier result suggests that the `-n` option of Sickle only removed 1,086 sequences. Can double check that `-n` option removed all sequences with at least one N:

```bash
awk 'NR % 4 == 2' 2_10_AGTTCC_L008_R1_001_trimmedx2_no_n.fastq  | grep -c N
0
```

Not quite sure what is happening here, and in any case we probably don't want to use `-n` option because most occurences of an N are just where there is 1 N in the read. These may still be useful for aligning to a reference.

## Testing Scythe + Sickle combined

With some difficulty, I made a make file that runs Scythe and Sickle on a small number of test FASTQ files. This makefile (`/share/tamu/Code/run_scythe_and_sickle.mk`) removes intermediate output files, leaving only processed FASTQ files.

The details of testing will not be included here.

# Main run of Sickle and Scythe #

Make new directory which will contain symbolic links to all FASTQ files and then analysis output files will be created in this directory from subsequent steps.

```bash
cd /share/tamu/Analysis
mkdir All_FASTQ_files
cd All_FASTQ_files
find /share/tamu/Data/DNase-Seq/Mouse/ -type f -name "*.fastq.gz" -exec ln -s {} \;
```

I made a simple bash script (`/share/tamu/Code/run_scythe_and_sickle.sh`) which is really just a wrapper around the make file. Can then submit this job to the queue with qsub command:

```bash
qsub -S /bin/bash -pe threaded 2 -M keith@bradnam.co /share/tamu/Code/run_scythe_and_sickle.sh
```

This ran overnight and finished creating a *_processed.fastq file for each of the 120 input files. Two of the output files were empty (`8Lvr_20_GTCCGC_L007_R1_001_processed.fastq` and `8Lvr_50_GTGAAA_L007_R1_001_processed.fastq`). Checking these by hand, I see that the input *fastq.gz files are also empty, so we can just ignore these from now on;

```bash
rm 8Lvr_20_GTCCGC_L007_R1_001_processed.fastq 8Lvr_50_GTGAAA_L007_R1_001_processed.fastq
```

# Testing short read mappers #

## Installing bowtie2 (version 2.2.4) ##

1. Cloned locally from https://github.com/BenLangmead/bowtie2.git
2. `mkdir -p Packages/Bowtie2`
3. Use `scp` to copy local `scythe` directory to `Packages/Bowtie2/1d0bc49125` # use SHA-1 hash for directory name
4. `cd Packages/Bowtie2`
5. `ln -s 1d0bc49125 current`
6. `cd current`
8. `make all`
8. No errors reported at this point, so…
9. `cd /share/tamu/bin`
10. `ln -s ../Packages/Bowtie2/current/bowtie2`
11. `ln -s ../Packages/Bowtie2/current/bowtie2-build`


## Installing BWA (version 0.7.10) ##

1. Downloaded from http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download
2. `mkdir -p Packages/BWA`
3. Use `scp` to copy downloaded BWA directory to `Packages/BWA/`
4. `cd Packages/BWA`
5. `bunzip2 bwa-0.7.10.tar.bz2`
6. `tar -xvf bwa-0.7.10.tar`
7. `ln -s bwa-0.7.10 current` 
7. `cd current`
8. `make`
8. No errors reported at this point, so…
9. `cd /share/tamu/bin`
10. `ln -s ../Packages/BWA/current/bwa`

## Installing Aryana (version 0.1) ##

1. Cloned locally from https://github.com/aryana-aligner/src.git
2. `mkdir -p Packages/Aryana`
3. Use `scp` to copy local `scythe` directory to `Packages/Aryana/a2712c0e6b` # use SHA-1 hash for directory name
4. `cd Packages/Aryana`
5. `ln -s a2712c0e6b current`
6. `cd current`
8. `make`
8. No errors reported at this point, so…
9. `cd /share/tamu/bin`
10. `ln -s ../Packages/Aryana/current/aryana`


## Testing Bowtie 2 ##

Make test directory and set up links to some reads and target genome (will just use a really small chromosome file for testing): 

```bash
mkdir /share/tamu/Analysis/Test/Bowtie2_test
cd /share/tamu/Analysis/Test/Bowtie2_test
ln -s /share/tamu/Analysis/All_FASTQ_files/1b_10_TGACCA_L008_R1_001_processed.fastq
cp cp /share/tamu/Data/Genomes/Mouse/mm9/chromosomes/chr19.fa.gz .
gunzip chr19.fa.gz
```

First step is to use `bowtie2-build` to make an index from our reference sequence(s):
```bash
bowtie2-build chr19.fa test_index
```

For a small reference file, this took ~2.5 minutes and ends up producing:

```bash
s -l
total 150013
lrwxrwxrwx  1 keith keith       77 Nov 23 06:13 1b_10_TGACCA_L008_R1_001_processed.fastq -> /share/tamu/Analysis/All_FASTQ_files/1b_10_TGACCA_L008_R1_001_processed.fastq
-rw-r--r--+ 1 keith keith 62569286 Nov 23 07:33 chr19.fa
-rw-r--r--+ 1 keith keith 23575307 Nov 23 07:35 test_index.1.bt2
-rw-r--r--+ 1 keith keith 14535564 Nov 23 07:35 test_index.2.bt2
-rw-r--r--+ 1 keith keith       44 Nov 23 07:34 test_index.3.bt2
-rw-r--r--+ 1 keith keith 14535558 Nov 23 07:34 test_index.4.bt2
-rw-r--r--+ 1 keith keith 23575307 Nov 23 07:36 test_index.rev.1.bt2
-rw-r--r--+ 1 keith keith 14535564 Nov 23 07:36 test_index.rev.2.bt2
```

Main run options use -x for index prefix name, -U for filename(s) of unpaired read data, -p for number of processors, and -S for SAM output file name. So for our first test run:

```bash
time bowtie2 -x test_index -U 1b_10_TGACCA_L008_R1_001_processed.fastq -p 4 -S test_output.sam
3957919 reads; of these:
  3957919 (100.00%) were unpaired; of these:
    3367827 (85.09%) aligned 0 times
    126726 (3.20%) aligned exactly 1 time
    463366 (11.71%) aligned >1 times
14.91% overall alignment rate

real    1m48.851s
user    6m57.508s
sys     0m7.825s
``` 

By default, bowtie looks for multiple distinct, valid alignments but reports only the best one. Can use -k <N> option to produce results for N best matches or -a mode to report on all matches. Let's investigate difference in final output when we use -k = 2, 5, or 10 and the -a option:

```bash
time bowtie2 -x test_index -U 1b_10_TGACCA_L008_R1_001_processed.fastq -p 4 -k 2 -S test_output_k2.sam
3957919 reads; of these:
  3957919 (100.00%) were unpaired; of these:
    3365308 (85.03%) aligned 0 times
    126304 (3.19%) aligned exactly 1 time
    466307 (11.78%) aligned >1 times
14.97% overall alignment rate

real    1m24.787s


time bowtie2 -x test_index -U 1b_10_TGACCA_L008_R1_001_processed.fastq -p 4 -k 5 -S test_output_k2.sam
3957919 reads; of these:
  3957919 (100.00%) were unpaired; of these:
    3362694 (84.96%) aligned 0 times
    126127 (3.19%) aligned exactly 1 time
    469098 (11.85%) aligned >1 times
15.04% overall alignment rate

real    1m51.474s


time bowtie2 -x test_index -U 1b_10_TGACCA_L008_R1_001_processed.fastq -p 4 -k 10 -S test_output_k10.sam
3957919 reads; of these:
  3957919 (100.00%) were unpaired; of these:
    3361257 (84.92%) aligned 0 times
    126034 (3.18%) aligned exactly 1 time
    470628 (11.89%) aligned >1 times
15.08% overall alignment rate

real    2m28.565s

```

I abandoned the -a option after about 20 hours as it hadn't finished and the output SAM file had grown to over 12 GB. Increasing value of -k, also resuts in larger output file sizes:

```bash
ls -lh *.sam
-rw-r--r--+ 1 keith keith 683M Nov 23 07:41 test_output.sam
-rw-r--r--+ 1 keith keith 794M Nov 23 07:57 test_output_k2.sam
-rw-r--r--+ 1 keith keith 1.1G Nov 23 08:00 test_output_k5.sam
-rw-r--r--+ 1 keith keith 1.6G Nov 23 08:02 test_output_k10.sam
-rw-r--r--+ 1 keith keith  12G Nov 24 04:04 test_output_a.sam
```



## Make Bowtie 2 index for all Mouse chromosome files ##

```bash
qlogin 
cd /share/tamu/Data/
mkdir Bowtie2_indexes
cd Bowtie2_indexes/
cp /share/tamu/Data/Genomes/Mouse/mm9/chromosomes/*.fa.gz .
gunzip *.fa.gz
time bowtie2-build chr10.fa,chr11.fa,chr12.fa,chr13.fa,chr13_random.fa,chr14.fa,chr15.fa,chr16.fa,chr16_random.fa,chr17.fa,chr17_random.fa,chr18.fa,chr19.fa,chr1.fa,chr1_random.fa,chr2.fa,chr3.fa,chr3_random.fa,chr4.fa,chr4_random.fa,chr5.fa,chr5_random.fa,chr6.fa,chr7.fa,chr7_random.fa,chr8.fa,chr8_random.fa,chr9.fa,chr9_random.fa,chrM.fa,chrUn_random.fa,chrX.fa,chrX_random.fa,chrY.fa,chrY_random.fa mm9_index

rm -f *.fa
cp /share/tamu/Data/Genomes/Mouse/mm10/chromosomes/*.fa.gz .
gunzip *.fa.gz
cat *.fa > all.fa
time bowtie2-build all.fa mm10_index
rm -f *.fa
exit
```

Takes about 4 hours for each version of the genome


## Test of Bowtie against full (mm9 and mm10) genomes ##

If search for 1 FASTQ file against 1 chromosome only takes ~2 minutes, now want to see how this changes against the full genome:

```bash
bowtie2 -x /share/tamu/Data/Bowtie2_indexes/mm9_index -U 1b_10_TGACCA_L008_R1_001_processed.fastq -p 4 -S test_output_mm9.sam
3957919 reads; of these:
  3957919 (100.00%) were unpaired; of these:
    182533 (4.61%) aligned 0 times
    2694414 (68.08%) aligned exactly 1 time
    1080972 (27.31%) aligned >1 times
95.39% overall alignment rate

bowtie2 -x /share/tamu/Data/Bowtie2_indexes/mm10_index -U 1b_10_TGACCA_L008_R1_001_processed.fastq -p 4 -S test_output_mm10.sam
3957919 reads; of these:
  3957919 (100.00%) were unpaired; of these:
    174969 (4.42%) aligned 0 times
    2694299 (68.07%) aligned exactly 1 time
    1088651 (27.51%) aligned >1 times
95.58% overall alignment rate
```


Surprisingly, this only took about 7.5 minutes for the 1 FASTQ file. Note that mm10 has slightly higher overall alignment rate.

Now let's see what difference the `--very-sensitive` option makes (this is a synonym for `-D 20 -R 3 -N 0 -L 20 -i S,1,0.50` options) and should be slower.

```bash
bowtie2 -x /share/tamu/Data/Bowtie2_indexes/mm9_index -U 1b_10_TGACCA_L008_R1_001_processed.fastq -p 4 --very-sensitive -S test_output_mm9_vs.sam
3957919 reads; of these:
  3957919 (100.00%) were unpaired; of these:
    161035 (4.07%) aligned 0 times
    2670520 (67.47%) aligned exactly 1 time
    1126364 (28.46%) aligned >1 times
95.93% overall alignment rate

bowtie2 -x /share/tamu/Data/Bowtie2_indexes/mm10_index -U 1b_10_TGACCA_L008_R1_001_processed.fastq -p 4 --very-sensitive -S test_output_mm10_vs.sam
3957919 reads; of these:
  3957919 (100.00%) were unpaired; of these:
    155621 (3.93%) aligned 0 times
    2670058 (67.46%) aligned exactly 1 time
    1132240 (28.61%) aligned >1 times
96.07% overall alignment rate
```

This takes about 11 minutes per run (compared to 7.5 for default parameters) and increases the overall alignment rate (95.39 -> 95.93% for mm9, 95.58 -> 96.07% for mm10).



### Distribution of mapping quality scores in mm10 ###

SAM format uses scores from 0 to 42 for mapping quality, where the score is a function of the probability that the mapping was wrong (though not all tools may be calculating probability in a way that makes sense, see [this post](http://biofinysics.blogspot.com/2014/05/the-slow-death-of-term-uniquely.html). Higher scores are better.

```bash
test_output_mm10.sam | cut -f 5 | sort | uniq -c | sort -nk 2
 260715 0
 574893 1
   8584 2
   8729 3
   1789 4
   1927 5
  14425 6
  18143 7
  10330 8
   2002 11
    440 12
   3182 14
   9048 15
   5756 16
    132 17
    243 18
    336 21
   4662 22
   8339 23
  37122 24
    349 25
   4950 26
   2251 27
 140363 30
   3971 31
  24658 32
    199 33
  10203 34
  26707 35
  19488 36
  33489 37
  33667 38
  63795 39
  15758 40
2607274 42
```

Not sure if we will want to filter out some of the matches with scores < 42.
From a [great blog post](http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html) on filtering MAPQ scores in Bowtie 2 output:

>If someone wanted to exclude "true multireads" from their data set, using MAPQ >= 2 would work. This would also exclude any uniquely mapping reads with >=4 mismatches over high quality bases. In terms of high quality bases and unireads, MAPQ >= 3 allows up to 3 mismatches, MAPQ >= 23 allows up to 2 mismatches, MAPQ >= 40 allows up to 1 mismatch, and MAPQ >= 42 allows 0 mismatches. There will also be other "maxireads" in most or all of these sets.


### Bowtie output options to consider using ###

`--no-unal` suppress reads that don't align
`--no-hd` suppress SAM header lines
`--omit-sec-seq` replace sequence and quality fields with asterisks (but only for secondary alignments)


# Main run of Bowtie 2 against mm9 and mm10 #

I will use a Makefile strategy again which will generate 2 SAM files for each input FASTQ file (1 for mm9 and 1 for mm10). Can't decide whether to keep information on non-aligning reads (as a) I won't be using them for anything and b) they take up disk space) or just omit. Will still need some post-processing steps to decide what to do with reads with low MAPQ scores.

I made a simple bash script (`/share/tamu/Code/run_bowtie2.sh`) which is just a wrapper around the `run_bowtie2.mk` make file. Can then submit this job to the queue with qsub command:

```bash
qsub -S /bin/bash -pe threaded 2 -M keith@bradnam.co /share/tamu/Code/run_bowtie2.sh
```

This didn't seem to work as planned. In the end I made a simple Perl wrapper script (`/share/tamu/Code/run_bowtie2.pl`) to submit all 236 jobs as single-threaded bowtie2 runs. This filled up the cluster queue with 1 job per file to be processed. It seemed to start very slowly, but sped up in the end (finished within ~2 days). As well as the main SAM file for each processed FASTQ file, this script also captures standard error from Bowtie 2 into *.err files (these contain summary mapping statistics). There is also an output file from the qsub process (`bowtie2.*`) which were all empty as output files were specified as part of the `bowtie2` command.

## Cleaning up ##

After this step, I could gzip the processed fastq files from earlier, and remove the empty bowtie2 qsub output files.

```bash
rm -f bowtie2.*
gzip *processed.fastq
```


# Queue management tips #

Request a node with 32G of memory: `qlogin -l h_vmem=32G`

See your jobs that are running: `qstat`
See your jobs plus information for all nodes: `qstat -f`
As above but show details for all users: `qstat -f -u '*'`
See summary of cluster queues: `qstat -g c`

List parallel environments that are available: `qconf -spl`


# Generating summary output from bowtie2 #

Want a Perl script that will:

+ process all *.err files and combine read/mapping counts for each sample (represented by several output files)
+ process SAM files to maybe make 2–3 counts of mapped reads based on MAPQ score (with sucessive levels of filtering)

Want to know whether any sample has a particular bias:

+ in total number of reads
+ in percentage of reads mapped

Also want to know how different mapping results are to mm9 vs mm10.

```bash
qlogin
generate_bowtie2_summary_stats.pl > bowtie_summary_stats.tsv
```

# Quality control #

A range of MAPQ scores are present (see earlier section) between zero (worst score possible) and 42 (best score possible). From inspection of SAM file (particularly alignment score field (AS:i), alternative alignment score (XS:i), and number of mismatches (XM:i), you can start to get a feel for what these lower MAPQ scores relate to.

+ You can have 1-3 mismatches but still get a MAPQ score of 42 (alignment score drops from zero to -2 through to -6)
+ MAPQ scores at 39 or lower start to have secondary alignments (XS:i field is present)
+ MAPQ scores >23 never have more than 3 mismatches
+ MAPQ scores of 23 can have 1–5 mismatches
+ MAPQ scores of 2 or lower can have 6-9 mismatches

Will keep alignments of 24 of better and will also combine SAM output from multiple SAM files from each sample. Finally, we will just use mm10 for now, rather than duplicating everything. Use a Perl script:

```bash
/share/tamu/Code/filter_sam_files.pl 24
Discarded 82308088 out of 416642306 alignments with MAPQ scores < 24
```

# Tidy up #

Can now move the processed SAM files into their own directory and zip the original SAM files.

```bash
cd /share/tamu/Analysis
mkdir SAM_files
cd SAM_files
mv ../All_FASTQ_files/*MAPQ24.sam .
cd ../share/tamu/Analysis/All_FASTQ_files
gzip *.sam
```