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

### Testing Scythe: part 1 ###
1. Using supplied sample files that are part of Scythe installation
2. `mkdir -p Analysis/Scythe_test`
3. `cd Analysis/Scythe_test`
4. Run Scythe:

```bash
scythe -a ../../Packages/Scythe/current/illumina_adapters.fa -o trimmed_seqs.fastq ../../Packages/Scythe/current/testing/reads.fastq
prior: 0.300

Adapter Trimming Complete
contaminated: 1329, uncontaminated: 8671, total: 10000
contamination rate: 0.132900
```

Resulting output file does indeed to be trimmed:

```bash
ls -lh trimmed_seqs.fasta ../../Packages/Scythe/current/testing/reads.fastq
-rw-r--r--+ 1 keith dglemay 2.3M Oct 17 06:45 ../../Packages/Scythe/current/testing/reads.fastq
-rw-r--r--+ 1 keith dglemay 2.0M Oct 17 07:02 trimmed_seqs.fastq
```


### Extracting adaptor sequences ###

From supplied Illumina PDF file, I made a single FASTA file containing sequence of the Universal TruSeq adapter, and the 24 TruSeq index adapters, and placed in  `Data/DNase-Seq/Mouse/adapter_sequences.fasta`

I then extracted the barcodes from the index adapters as follows:

```bash
cat adapter_sequences.fasta | sed 's/GATCGGAAGAGCACACGTCTGAACTCCAGTCAC//' | grep -v ">" | tail -n 24 | cut -c 1-6 > barcodes.txt
```

Presumably barcodes in file names should match barcodes in the FASTQ identifiers in the respective files. But this is not what happens:

```bash
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


### Testing Scythe: part 2 ###

Now to try with some real DNAse-Seq data (which is still all zipped up).

1. `cd Analysis/Scythe_test`
2. `cp ../../Data/DNase-Seq/Mouse/Sample_2_10.zip .` - this is the smallest file
3. `unzip Sample_2_10.zip`
4. `cd Sample_2_10`

Main run:
```bash
scythe -a ../../../Data/DNase-Seq/Mouse/adapter_sequences.fasta -o 2_10_AGTTCC_L008_R1_001_trimmed.fastq 2_10_AGTTCC_L008_R1_001.fastq.gz

prior: 0.300

Adapter Trimming Complete
contaminated: 264003, uncontaminated: 3735997, total: 4000000
contamination rate: 0.066001

```