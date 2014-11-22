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

## Installing Aryana (version ?) ##

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