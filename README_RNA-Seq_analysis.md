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

