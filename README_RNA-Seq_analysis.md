Milk-DNase-Seq-Project: RNA-Seq Analyis
=======================================

See main [README](README.md) file for more information about this project.

## Bovine RNA-seq data ##

Stored in `/share/tamu/Data/RNA-Seq/Cow/2014-10`. Looks like paired-read 100 bp data. In total 31 x 2 files, ranging from 1–3.5 GB in size. See also the `/share/tamu/Data/RNA-Seq/Cow/Metadata` directory which contains a metadata file which suggests that we have data from 15 virgin cows and 16 lacating cows.

The ultimate goal is to find genes that are differentially expressed between these two developmental stages.

These files were originally compressed with bzip2, will re-compress with gzip so that existing pipelines can work with them. And will also rename them to have fastq suffix:

```bash
cd /share/tamu/Data/RNA-Seq/Cow/2014-10
bunzip2 *.bz2
rename.pl s/txt/fastq/ *.txt
gzip *.fastq
```



## Checking barcodes in RNA-Seq data ##

Let's check on all barcodes being used. Will make some soft links to the data files

```bash
cd Analysis/Test
mkdir RNA-Seq_Barcode_check
cd RNA-Seq_Barcode_check

qlogin

bunzip2 -c ../../../Data/RNA-Seq/Cow/2014-10/*.bz2 | grep "@HWI" | sed 's/.* [12]:N:0://' | sort | uniq -c > barcodes_in_identifiers.txt
```

Unfortunately this failed due to a 'No space left on device error'. So maybe need to treat each file separately.


# Test run of Scythe and Sickle#

Unlike the DNase-Seq data, we now have paired-end data, which requires running Sickle a little differently. So first, let's do a test (using 10,000 reads from each of two paired FASTQ files):

```
cd /share/tamu/Analysis/Test
mkdir Paired_end_scythe_sickle_test
cd Paired_end_scythe_sickle_test

gunzip -c /share/tamu/Data/RNA-Seq/Cow/2014-10/L468_C4K66ACXX-7-ID09_1_sequence.fastq.gz | head -n 400000 > test_100K_1.fastq

gunzip -c /share/tamu/Data/RNA-Seq/Cow/2014-10/L468_C4K66ACXX-7-ID09_2_sequence.fastq.gz | head -n 400000 > test_100K_2.fastq

scythe -a /share/tamu/Data/DNase-Seq/Mouse/adapter_sequences.fasta -o test_100K_1_scythe.fastq test_100K_1.fastq

scythe -a /share/tamu/Data/DNase-Seq/Mouse/adapter_sequences.fasta -o test_100K_2_scythe.fastq test_100K_2.fastq
```

Scythe does not need to do anything for paired-end files, just run against each file separately. This gives output like so:

```
Adapter Trimming Complete
contaminated: 7131, uncontaminated: 92869, total: 100000
contamination rate: 0.071310
```

Sickle has a paired end mode which is run slightly differently than with unpaired data (have up to 3 output files):

```
sickle pe -f test_100K_1_scythe.fastq -r test_100K_2_scythe.fastq -t sanger -o test_100K_1_sickle.fastq -p test_100K_2_sickle.fastq -s test_100K_single_sickle.fastq

PE forwrd file: test_100K_1_scythe.fastq
PE reverse file: test_100K_2_scythe.fastq

Total input FastQ records: 200000 (100000 pairs)

FastQ paired records kept: 194482 (97241 pairs)
FastQ single records kept: 2700 (from PE1: 2592, from PE2: 108)
FastQ paired records discarded: 118 (59 pairs)
FastQ single records discarded: 2700 (from PE1: 108, from PE2: 2592)

wc -l *.fastq
   400000 test_100K_1.fastq
   400000 test_100K_1_scythe.fastq
   388964 test_100K_1_sickle.fastq
   400000 test_100K_2.fastq
   400000 test_100K_2_scythe.fastq
   388964 test_100K_2_sickle.fastq
    10800 test_100K_single_sickle.fastq  
```

But is this different from running sickle in single read mode (se)?

```bash
sickle se -f test_100K_1_scythe.fastq -t sanger -o test_100K_1_sickle_se.fastq

SE input file: test_100K_1_scythe.fastq

Total FastQ records: 100000
FastQ records kept: 99833
FastQ records discarded: 167

sickle se -f test_100K_2_scythe.fastq -t sanger -o test_100K_2_sickle_se.fastq

SE input file: test_100K_2_scythe.fastq

Total FastQ records: 100000
FastQ records kept: 97349
FastQ records discarded: 2651

wc -l *sickle*fastq
   388964 test_100K_1_sickle.fastq
   399332 test_100K_1_sickle_se.fastq
   388964 test_100K_2_sickle.fastq
   389396 test_100K_2_sickle_se.fastq
    10800 test_100K_single_sickle.fastq
```

Seemingly, you can't run sickle in single read mode on paired files without ending up with different numbers of reads in each file of the pair (which could be problematic for downstream steps that need paired data).

# Main run of Scythe and Sickle#

Make new directory which will contain symbolic links to all FASTQ files and then analysis output files will be created in this directory from subsequent steps.

```bash
cd /share/tamu/Analysis
mkdir RNA-Seq_FASTQ_files
cd RNA-Seq_FASTQ_files
find /share/tamu/Data/RNA-Seq/Cow/2014-10/ -type f -name "*.fastq.gz" -exec ln -s {} \;
```

I made a simple bash script (`/share/tamu/Code/run_scythe_and_sickle2.sh`) which simply runs scythe on each FASTQ file from a pair, and then runs sickle (in paired end mode) on both output files of scythe. 

Then I used a Perl script (`wrapper2.pl`) to automate the submission of all of these jobs with qsub. Thirty-one jobs in total to process 62 paired FASTQ files.


# Installing TopHat

TopHat is [available here](http://ccb.jhu.edu/software/tophat/index.shtml).

```bash
cd /share/tamu/Packages
mkdir TopHat
cd TopHat
curl -O http://ccb.jhu.edu/software/tophat/downloads/tophat-2.0.13.Linux_x86_64.tar.gz
tar -xvzf tophat-2.0.13.Linux_x86_64.tar.gz
ln -s tophat-2.0.13.Linux_x86_64 current
cd /share/tamu/bin
ln -s ../Packages/TopHat/current/tophat
ln -s ../Packages/TopHat/current/tophat2
```

# Make Bowtie2 index

Will make indexes for UMD 3.1 (bosTau6) and newer BCM Btau_4.6.1 (bosTau7). But start with the genome that we will probably use the most. We also want to have same FASTA file present in the index directory (for use by TopHat later on). Will try just using symbolic link for now:

```bash
cd /share/tamu/Data/Bowtie2_indexes
time bowtie2-build ../Genomes/Cow/bosTau6/bigZips/bosTau6.fa bosTau6_index
time bowtie2-build ../Genomes/Cow/bosTau6/Ensembl-78-genome/Bos_taurus.UMD3.1.dna.toplevel.fa Ensembl-78-index
ln -s ../Genomes/Cow/bosTau6/bigZips/bosTau6.fa bosTau6_index.fa
ln -s ../Genomes/Cow/bosTau7/bigZips/bosTau7.fa bosTau7_index.fa
ln -s ../Genomes/Cow/bosTau6/Ensembl-78-genome/Bos_taurus.UMD3.1.dna.toplevel.fa Ensembl-78-index.fa
```



# Main run on TopHat

After [extensive testing](README_TopHat_testing.md) I chose the following parameters for TopHat runs:

```bash
# run39: standard deviation = 1000, -r = 800 and use --b2-very-sensitive option
time tophat --no-convert-bam -o run39 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --b2-very-sensitive --no-mixed --no-discordant -T -M -r 800 --mate-std-dev 1000 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq
```

Use `wrapper3.pl` script in `/share/tamu/Analysis/RNA-Seq_FASTQ_files` directory to submit `run_tophat.sh` script to job scheduler. This should create new output directories in `/share/tamu/Analysis/TopHat_output` (one directory for each of the 31 samples). 

Reminder, this run is not going to use the singles files. We may wish to return to these later when looking to improve gene annotations.

All but one job failed with 'MemoryError' so I resubmitted with a request for 8 GB RAM per thread and instructed TopHat and qsub to just use one thread per process. Qsub commands looked like this:

```bash
qsub -S /bin/bash -pe threaded 1 -l h_vmem=8G -M keith\@bradnam.co -m be -N krb_tophat_run /share/tamu/Code/run_tophat.sh L478_C4K66ACXX-7-ID10_1_sequence_processed.fastq L478_C4K66ACXX-7-ID10_2_sequence_processed.fastq /share/tamu/Analysis/TopHat_output/L478_C4K66ACXX-7-ID10
```


# Generate count data

Now we have SAM output files, we need to generate counts of reads mapped to each transcript. We also need to do some filtering of the SAM files:

1. Only want uniquely mapped reads
2. Only want concordantly mapped pairs of reads
3. Only want read pairs that map to the same chromosome
4. Might only want to keep reads that not too far apart (< 1 Mbp?)

## Install HTSeq

We will use [HTSeq](http://www-huber.embl.de/users/anders/HTSeq/) software (a Python package) for generating read counts. Need to install this on isner.

```bash
cd /Chanlab/Packages
mkdir HTSeq
cd HTSeq
curl -O https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1.tar.gz
tar -xvzf HTSeq-0.6.1.tar.gz
ln -s current HTSeq-0.6.1
cd current
python setup.py build
sudo python setup.py install
cd /Chanlab/bin
ln -s ../Packages/HTSeq/current/scripts/htseq-count
ln -s ../Packages/HTSeq/current/scripts/htseq-qa
```

## Filtering SAM files ##

Start off by making links to pre-existing SAM files as well as making a couple of test files (1 or 5 million reads)

```bash
cd /Chanlab/Scratch/keith/RNAseq
mkdir Filter_SAM_files
cd Filter_SAM_files
grep -v "^@"  ../TopHat_output/L468_C4K66ACXX-7-ID09/accepted_hits.sam | head -n 1000000 > test_1M.sam
grep -v "^@"  ../TopHat_output/L468_C4K66ACXX-7-ID09/accepted_hits.sam head -n 5000000 > test_5M.sam
grep -v "^@"  ../TopHat_output/L468_C4K66ACXX-7-ID09/accepted_hits.sam head -n 25000000 > test_25M.sam

```

Two ways to get only the good quality reads (unique and concordantly mapped). Both of the following give identical output, but the awk method doesn't require SAM header lines to be present:

```bash
cat test_1M.sam | awk '$5 == 50 && ($2 == 163 || $2 == 147 || $2 == 83 || $2 == 99) {print}' 

samtools  view -q 50 -f 0x2 test_1M.sam


cat test_1M.sam | awk '$5 == 50 && ($2 == 163 || $2 == 147 || $2 == 83 || $2 == 99) {print}' > test_1M_HQ.sam

cat test_5M.sam | awk '$5 == 50 && ($2 == 163 || $2 == 147 || $2 == 83 || $2 == 99) {print}' > test_5M_HQ.sam

cat test_25M.sam | awk '$5 == 50 && ($2 == 163 || $2 == 147 || $2 == 83 || $2 == 99) {print}' > test_25M_HQ.sam

```

# Testing HTSeq

Need to use `-r pos` to tell htseq-count that our data is already sorted by position. The default GFF feature to map to is exon (`-t` option) which should suffice as should the default value of 'gene_id' for the `-i` option (specifying the attribute field). First runs will test difference of using `-s no` to specify that this is not a stranded assay (default is yes).

```bash

htseq-count -r pos -m intersection-nonempty test_1M_HQ.sam ../Bowtie2_indexes/Ensembl_78_transcriptome.gff > test_htseq_run1.txt

htseq-count -r pos -m intersection-nonempty -s no test_1M_HQ.sam ../Bowtie2_indexes/Ensembl_78_transcriptome.gff > test_htseq_run2.txt

```

We definitely need to use `-s no` otherwise it doesn't count reads that are on the opposite strand to an exon feature.

# Running HTseq

Wrote a simple wrapper Perl script to first filter SAM files (using awk approach shown above) and then run htseq-count. Final output files are in `/share/tamu/Analysis/Filter_SAM_files`

I then wrote another Perl wrapper script `/share/tamu/Code/join_htseq_count_files.pl` that uses the Unix `join` command to join all of the counts files into one final file:

```bash
join_htseq_count_files.pl "*-ID*.txt"
mv tmp_join_file.txt count_data_for_all_samples.tsv
```





