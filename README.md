Milk-DNase-Seq-Project
======================

Analysis of milk DNase-Seq data in mouse and cow.


## Setup of directories for raw data

In the shared data folder, we will create:

	Data/
	Data/Genomes/
	Data/Genomes/Cow/
	Data/Genomes/Mouse/	
	Data/Genomes/Mouse/mm9/
	Data/Genomes/Mouse/mm10/
	Data/DNAse-Seq/
	Data/DNAse-Seq/Cow/
	Data/DNAse-Seq/Mouse/
	Data/RNA-Seq/
	Data/RNA-Seq/Cow/
	Data/RNA-Seq/Mouse/

	
`Genomes` subdirectory will contain reference genome sequences and annotations, whereas
`DNAse-Seq` and `RNA-Seq` directories will contain our raw experimental data (in date-
versioned subdirectories).
 
I've also added a `Data/RNA-Seq/Cow/Metadata` for a couple of files Danielle provided
about the Cow RNA-Seq data:

`CowLactationRNASeqMetadata.txt` - A guess of what the sample IDs really mean
`CowLactationRNASeqMappingStats.tsv` - tab delimited file containing mapping statistics from Baylor