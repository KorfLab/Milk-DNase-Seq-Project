Milk-DNase-Seq-Project
======================

Analysis of milk DNase-Seq data in mouse and cow.


## Setup of directories for raw data ##

In the shared folder, we will create:

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

### Testing Scythe ###
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
