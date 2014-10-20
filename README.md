Milk-DNase-Seq-Project
======================

Analysis of milk DNase-Seq data in mouse and cow.


## Setup of directories for raw data ##

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

### Genome data ###
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

### Extracting adaptor sequences ###

From supplied Illumina PDF file, I made a single FASTA file containing sequence of the Universal TruSeq adapter, and the 24 TruSeq index adapters, and placed in  `Data/DNase-Seq/Mouse/adapter_sequences.fasta`

I then extracted the barcodes from the index adapters as follows:

```bash
cat adapter_sequences.fasta | sed 's/GATCGGAAGAGCACACGTCTGAACTCCAGTCAC//' | grep -v ">" | tail -n 24 | cut -c 1-6 > barcodes.txt
```


## Other directories for project ##

A top level `Code` folder will contain scripts used to analyze data. This may be expanded to contain subdirectories (if we have some code in Git repositories).

A top level `Packages` folder will be used to install pre-existing bioinformatics tools. Until `git` is available on the cluster, may have to `git clone` locally and then `scp` files to our machine on the cluster.

A top level `bin` directory will contain symbolic links to executables in either `Code` or `Packages`. This `bin` directory should be added to users' $PATH.

A top level `Analysis` directory which will contain results from any steps we run.

## Installing Scythe (v0.991) ##
1. Cloned locally from https://github.com/vsbuffalo/scythe
2. `mkdir -p Packages/Scythe`
3. Use `scp` to copy local `scythe` directory to `Packages/Scythe/c128b19c65` # use SHA-1 hash for directory name
4. `cd Packages/Scythe`
5. `ln -s c128b19c65 current`
6. `cd current`
7. `make all`
8. No errors reported at this point, so…
9. `cd /Share/tamu/bin`
10. `ln -s ../Packages/Scythe/current/scythe`

### Making test dataset ###

While testing Scythe and other tools, it will be good to just one work with a few FASTQ files. Make copy of `Sample_2_10.zip` sequences as this is the smallest file:

```bash
cd Analyis
mkdir Test
cp ../../Data/DNase-Seq/Mouse/Sample_2_10.zip .
unzip Sample_2_10.zip
cd Sample_2_10
```

### Checking barcodes in test dataset

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

### Testing Scythe: part 1 ###
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




### Testing Scythe: part 2 ###

Now to try with some real DNAse-Seq data (which is still all zipped up).

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