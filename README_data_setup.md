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
mkdir bosTau6 bosTau7 bosTau8
cd bosTau6
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/bosTau7/bigZips ./
cd bosTau6
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/bosTau7/bigZips ./
cd bosTau8
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/bosTau8/bigZips ./
```

Note that two versions of the mouse genome were downloaded and three versions of the cow genome (bosTau6–8). BosTau6 corresponds to UMD 3.1 (2009, with Ensembl annotations from 2011), BosTau7 corresponds to Btau_4.6.1 from BCM (2011?), and BosTau8 is the latest (2014) UMD 3.1.1 assembly. Also, for mm9 I downloaded the comparison alignments to bosTau7. Initially, just the 'chromosomes' subdirectories of the mouse URLs were copied.

Also note from UCSC Genome site:

>We are pleased to announce the release of a Genome Browser for the June 2014 assembly of cow, Bos taurus (Bos_taurus_UMD 3.1.1, UCSC version bosTau8). This updated cow assembly was provided by the UMD Center for Bioinformatics and Computational Biology (CBCB). This assembly is an update to the previous UMD 3.1 (bosTau6) assembly. UMD 3.1 contained 138 unlocalized contigs that were found to be contaminants. These have been suppressed in UMD 3.1.1.

Suggests that there may not be a huge difference in the utility of BosTau8 over BosTau6. I also copied the Ensembl version of what should be the same data:

```bash
cd /share/tamu/Data/Genomes/Cow/bosTau6
mkdir Ensembl-78-genome
cd Ensembl-78-genome/
curl -O ftp://ftp.ensembl.org//pub/release-78/fasta/bos_taurus/dna/Bos_taurus.UMD3.1.dna.toplevel.fa.gz
mkdir Ensembl-78-genome
cd Ensembl-78-genome/
curl -O ftp://ftp.ensembl.org//pub/release-78/fasta/bos_taurus/dna/README
```


## Cow genome annotation data ##

Ensembl is basing annotation on UMD 3.1 (BosTau6), so we should grab that info. Start with just the GTF data:

```bash
cd /share/tamu/Data/Genomes/Cow/bosTau6
mkdir Ensembl-78-annotations
cd Ensembl-78-annotations
curl -O ftp://ftp.ensembl.org//pub/release-78/gtf/bos_taurus/Bos_taurus.UMD3.1.78.gtf.gz
curl -O ftp://ftp.ensembl.org//pub/release-78/gtf/bos_taurus/README
gunzip Bos_taurus.UMD3.1.78.gtf.gz

```


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

