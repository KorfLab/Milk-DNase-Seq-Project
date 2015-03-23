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



# Basic test of TopHat #

Download the official [test dataset](http://ccb.jhu.edu/software/tophat/downloads/test_data.tar.gz)

```bash
mkdir /share/tamu/Analysis/Test/TopHat_test1
cd /share/tamu/Analysis/Test/TopHat_test1
curl -O http://ccb.jhu.edu/software/tophat/downloads/test_data.tar.gz
tar -xvzf test_data.tar.gz
mv test_data/* .
tophat -r 20 test_ref reads_1.fq reads_2.fq
```

## Test runs of TopHat to explore parameters ##

The above test with the simple test dataset seemed to work as expected. So can move to a bigger test with some of the real data (100K reads from each of paired end file, and 100K reads from the singles file):

```bash
mkdir /share/tamu/Analysis/Test/TopHat_test2
cd /share/tamu/Analysis/Test/TopHat_test2
head -n 400000 ../../RNA-Seq_FASTQ_files/L468_C4K66ACXX-7-ID09_1_sequence_processed.fastq > seqs_100K_1.fastq
head -n 400000 ../../RNA-Seq_FASTQ_files/L468_C4K66ACXX-7-ID09_2_sequence_processed.fastq > seqs_100K_2.fastq
head -n 400000 ../../RNA-Seq_FASTQ_files/L468_C4K66ACXX-7-ID09_sequence_singles.fastq > seqs_100K_singles.fastq
```

### Test run 1 ###
First step is to a minimal run without using singles file and not specifying any other options other than the minimum (index file and input read data). This will use a default value of 50 for the -r option that specifies distances between read pairs (I still don't know what this should be):

```bash
tophat --no-convert-bam /share/tamu/Data/Bowtie2_indexes/bosTau6_index seqs_100K_1.fastq seqs_100K_2.fastq
mv tophat_out run1
```

### Test run 2 ###

Test run 1 worked and produced output files in default output directory (subsequently renamed), now try adding in the singles file:

```bash
time tophat --no-convert-bam -o run2 /share/tamu/Data/Bowtie2_indexes/bosTau6_index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```

### Test run 3 ###

Now see if using Ensembl version of the genome gives same result as UCSC's version of the cow genome (should be the same sequence!)

```bash
time tophat --no-convert-bam -o run3 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```


### Test run 4 ###

Now try Ensembl version of the genome with annotations

```bash
time tophat --no-convert-bam -o run4b --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```

### Test run 5 ###

Now request concordant mappings only, but no Ensembl annotation

```bash
time tophat --no-convert-bam 
--no-discordant -o run5 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```


### Test run 6 ###

Now request --no-mixed option

```bash
time tophat --no-convert-bam 
--no-mixed -o run6 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```


### Test run 7 ###

Now request --no-mixed option and --no-discordant options

```bash
time tophat --no-convert-bam 
--no-mixed --no-discordant -o run7 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```


### Test run 8 ###

Does -r option make much difference (default insert size gap is 50), try 300 bp

```bash
time tophat --no-convert-bam 
--no-mixed --no-discordant -r 300 -o run8 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```

## Comparing output from runs ##

Adding in the singles sequences (run2) means you overall have many more reads that are mapped, though strangely it increased the number of mapped pairs that were mapped from run1. E.g. initially, 85371 of read pairs mapped in run1 , but the addition of the singles file increased this to 85704.

Run3 confirmed (thankfully!) that there is no different between the UCSC and Ensembl versions of the BosTau6 genome sequence. 

Run4 added the GTF file of annotations and this increased the overall mapping rate (89.7% in run3 to 92.5% in run4) as well as the concordant pair mapping rate (83.9% to 89.0%). Also many more reads mapped with max MAPQ score (271598 vs 254551).

The --no-discordant option (run5) results in slightly fewer aligned pairs mapped (83044 vs 85704 but produces slightly more mappings with max MAPQ score (255835 vs 254551). Also the align_summary file still implies that there are 0.9% discordant pairs vs 2.1% without the --no-discordant option.

The --no-mixed option (run6) ends up not mapping any of the single reads and when you add the --no-discordant option as well (run7) you still see a few discordant pairs remain, and it ends up mapping slightly fewer pairs.

Finally run8 used a different value for -r (insert pair distance) of 300 vs the default value of 50, but this was also still using the --no-mixed and --no-discordant options. This only resulted in *one* more pair being mapped, though the number of mapped read pairs with maximum MAPQ score increased (160286 to 160402).


## More testing ##

Want to focus on run4 (using GTF annotations) and explore use of other parameters relating to this:

### Test run 9 ###

Use -T option to only map to transcriptome from annotations (can compare results to run4 to see how many reads are matching outside the annotations):

```bash
time tophat --no-convert-bam -o run9 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf -T /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```


### Test run 10 ###

Use -M option to first pre-map to genome to get an idea of reads which might have multiple hits (genome-wide). From the TopHat manual:

>When mapping reads on the transcriptome, some repetitive or low complexity reads that would be discarded in the context of the genome may appear to align to the transcript sequences and thus may end up reported as mapped to those genes only. This option directs TopHat to first align the reads to the whole genome in order to determine and exclude such multi-mapped reads (according to the value of the -g/--max-multihits option).

```bash
time tophat --no-convert-bam -o run10 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf -M /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```


### Test run 11 ###

Use -M option to first pre-map to genome to get an idea of reads which might have multiple hits (genome-wide). But also use -T option to only map to transcriptome. Maybe this won't work, but can compare to run9 to see if -M option reduces number of aligned reads.

```bash
time tophat --no-convert-bam -o run11 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf -M -T /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```


### Test run 12 ###

One final test of -r option. This time want to keep other options simple. So effectively repeat run3 but just change -r to 200. Will use run3 and run12 data to assess length of inserts between correctly mapped read pairs.

```bash
time tophat --no-convert-bam -o run12 -r 200 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```

 

## Comparing output from runs (part 2) ## 

92.5% of reads mapped to transcriptome and genome (run4) whereas restricting matches to just the transcriptome (-T, run9) reduces this to 73.1%. Trying the pre-map option (-M) to first map to the genome before mapping to the genome (run10) made a slight increase to the mapping rate (up to 73.9%). When trying to combine -T and -M options (run11) we see similar mapping results to run9 (73.1%) but number of aligned read pairs is slightly lower (69679 vs 70800) as is concordant pair mapping rate (69.4% vs 70.6%).

The final run with an -r value of 200 to compare to the default value of 50 produced the same headline mapping rate (89.7%), the same number of aligned pairs (85,704) but a slightly higher rate of concordant mapped pairs (84.2% vs 83.9%). This suggests doing three more runs to test other values of -r:

### Test run 13 ###

-r = 100
```bash
time tophat --no-convert-bam -o run13 -r 100  /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```


### Test run 14 ###

-r = 300
```bash
time tophat --no-convert-bam -o run14 -r 300  /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```

### Test run 15 ###

-r = 400
```bash
time tophat --no-convert-bam -o run15 -r 400  /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```

### Effect of -r option ###

Number of concordantly mapped pairs by value of -r:

| r   | Mapping %    | Concordant pairs% | Concordant reads, MAPQ=50 | 
|-----|--------------|-------------------|---------------------------|
| 50  |    89.7%     |     83.9%         |         101,390           |
| 100 |    89.7%     |     84.0%         |         109,912           |
| 200 |    89.7%     |     84.2%         |         124,634           |
| 300 |    89.7%     |     84.0%         |         129,650           |
| 400 |    89.7%     |     84.3%         |         130,696           |



## Mapping quality values and discordant mapping ##

Inspecting MAPQ scores in resultant SAM files, about 85% of mapped reads (with default parameters) get MAPQ scores of 50 (unique matches). This means we may want to discard the ~15% of reads that have multiple matches (MAPQ scores of 0, 1, or 2). Also, if — for the purposes of calculating differential expression results — we only want to keep concordantly mapped pairs, we will want to keep results in SAM file that has a Bitwise SAM flag (column 2) of 99, 147, 163, or 83 (see [this post](https://ppotato.wordpress.com/2010/08/25/samtool-bitwise-flag-paired-reads/) for more info). This also means that we might not need to map singles file at this stage (but would return to this if we want to look at annotation issues). So this means we can run with --no-discordant option (and maybe --no-mixed option too).


## Double checking mapped insert sizes ##

Can use results from run3, run12, run13, run14, and run15 to estimate insert size from distant between concordantly mapped pairs (and work out whether changing -r option makes any notable difference to this). This can also help us calculate a standard deviation to use with --mate-std-dev option.


Want to see if we can identify insert size of mapped pairs from SAM data files. Want to ignore SAM header lines and only keep MAPQ scores of 50. For this purpose we may wish to also only focus on concordant mapped pairs (values 99, 147, 163, and 83 in column 2), and with perfect alignment scores (no negative AS: values). Because each read has a pair in the file, we can simply use `grep -v "-"` at one point (after cutting out certain columns first). With a bit of awk, we can get the actual insert size from the TLEN field (column 9), by subtracting twice the read length. However, as read lengths are not always 101 (imperfect mappings), we can cheat by only taking reads where all bases map (look for '101M' string in column 6). Altogether, this gives us something like:


```bash
cat run3/accepted_hits.sam | grep -vE "^@" | awk '$5 == 50 && ($2 == 99 || $2 == 147 || $2 == 163 || $2 == 83) {print}' | cut -f 6,9,12 | grep -v "-" | grep "101M" | cut -f 2 | awk '{print $1,"\t",$1-202}'  > run3_mapping_distances.tsv
```

Repeat for other runs in this series, and use R to make histograms from 2nd column of output.

### Test run 16 ###

As before in test run11, but just increase -r option to 400 to compare to run15. Will it make a difference to average insert sizes when you only map to the transcriptome?

```bash
time tophat --no-convert-bam -o run16 -r 400 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf -M -T /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```
 

### Test run 17 ###

As before in test run15, but use -max-std-dev option to see if this makes a difference (can use actual value we have from the data of run15 — 1382). The actual default value is 20!!!

```bash
time tophat --no-convert-bam -o run17 -r 400 --mate-std-dev 1382 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq,seqs_100K_singles.fastq
```


## Summary ##

Now we want to focus on using mappings to the transcriptome (-M -T option), but play around with different values of r and std deviation and focus on concordant mappings only so ignore singles file. Also will try data from a different individual.


### Test run 18 ###

Default value of -r and standard deviation. Map just to transcriptome (no -M option).

```bash
time tophat --no-convert-bam -o run18 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf --no-mixed --no-discordant -T /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq
```


### Test run 19 ###

As in run18 but add -M option 

```bash
time tophat --no-convert-bam -o run19 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf --no-mixed --no-discordant -T -M /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq
```


### Test run 20 ###

As in run19 but use -r 400

```bash
time tophat --no-convert-bam -o run20 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf --no-mixed --no-discordant -T -M -r 400 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq
```


### Test run 21 ###

As in run20 but use -r 800

```bash
time tophat --no-convert-bam -o run21 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf --no-mixed --no-discordant -T -M -r 800 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq
```


### Test run 22 ###

As in run21 but use --mate-std-dev of 1,000

```bash
time tophat --no-convert-bam -o run22 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf --no-mixed --no-discordant -T -M -r 800 --mate-std-dev 1000 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq
```

### Test run 23 ###

As in run22 but use --mate-std-dev of 2,000

```bash
time tophat --no-convert-bam -o run23 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf --no-mixed --no-discordant -T -M -r 800 --mate-std-dev 2000 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq
```

## Summary of recent test runs ##

Best = MAPQ score of 50, only proper concordant mapped pairs, and 101 bases
match (could still result in pairs of reads spanning splice junctions)

| Run | Notes       | -r  | Mapping% | Concordant% | MAPQ=50   | Best     |
|:---:|:-----------:|:---:|:--------:|:-----------:|:---------:|:--------:|
| 18  | T           |  50 |   68.9%  |    68.7%    |  136,748  |  84,996  |
| 19  | TM          |  50 |   69.0%  |    68.8%    |  136,730  |  40,071  |
| 20  | TM          | 400 |   69.0%  |    69.9%    |  137,026  |  53,109  |
| 21  | TM          | 800 |   69.0%  |    69.0%    |  136,910  |  54,545  |
| 22  | TM, std=1K  | 800 |   69.0%  |    69.0%    |  136,918  |  59,019  |
| 23  | TM, std=2K  | 800 |   69.0%  |    69.0%    |  136,916  |  60,019  |

Want to extract TLEN fields and convert to insert size, but maybe want to only take one read from a pair. Select for bitwise flags of 83 or 89 to get first reads?

```bash
cat run18/accepted_hits.sam | grep -vE "^@" | awk '$5 == 50 && ($2 == 99 || $2 == 83) {print}' |grep "101M" | cut -f 9 | grep '-' | awk '{print $1,"\t",-$1-202}'  > run18_mapping_distances.tsv
```

## Test run from different individuals ##

Maybe there is something very different with the test data from the selected individual. So try repeating run22 with data from 2 other cows:

```bash
head -n 400000 ../../RNA-Seq_FASTQ_files/L502_C5B2KACXX-7-ID02_1_sequence_processed.fastq > seqsB_100K_1.fastq
head -n 400000 ../../RNA-Seq_FASTQ_files/L502_C5B2KACXX-7-ID02_2_sequence_processed.fastq > seqsB_100K_2.fastq

head -n 400000 ../../RNA-Seq_FASTQ_files/L531_C5B2KACXX-7-ID06_1_sequence_processed.fastq > seqsC_100K_1.fastq
head -n 400000 ../../RNA-Seq_FASTQ_files/L531_C5B2KACXX-7-ID06_2_sequence_processed.fastq > seqsC_100K_2.fastq
```

### Test run 24 ###

As in run22 but with seqsB data

```bash
time tophat --no-convert-bam -o run24 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf --no-mixed --no-discordant -T -M -r 800 --mate-std-dev 2000 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqsB_100K_1.fastq seqsB_100K_2.fastq
```

### Test run 25 ###

As in run22 but with seqsC data

```bash
time tophat --no-convert-bam -o run25 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf --no-mixed --no-discordant -T -M -r 800 --mate-std-dev 2000 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqsC_100K_1.fastq seqsC_100K_2.fastq
```



## Bigger test run with final parameters

Run26: now try with 1 million sequences from each pair (10x previous runs, just under 1% of total file size of parent file: 100,101,576 reads).

```bash
head -n 4000000 ../../RNA-Seq_FASTQ_files/L468_C4K66ACXX-7-ID09_1_sequence_processed.fastq > seqs_1M_1.fastq
head -n 4000000 ../../RNA-Seq_FASTQ_files/L468_C4K66ACXX-7-ID09_2_sequence_processed.fastq > seqs_1M_2.fastq

time tophat --no-convert-bam -o run26 --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf --no-mixed --no-discordant -T -M -r 800 --mate-std-dev 1000 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_1M_1.fastq seqs_1M_2.fastq

```


## Build transcriptome index to save time in future ##

As a final check of inner distances, we want to map to transcriptome and keep coordinates relative to transcriptome only. Current mappings put things back in genome coordinates so read pairs that span introns have much larger inner distances. Can make TopHat generate the Transcriptome index files (that it otherwise builds and then deletes) and this can be used for all subsequent runs which will speed things up. Need to use the `--transcriptome-index <dir/prefix>` option which will first build index in specified directory. Subsequent runs of TopHat can then use the same option to reuse the index, and this means you no longer need the --GTF option. This will also let us calculate a more accurate value for --mate-std-dev.


```bash
time tophat --GTF /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index
```


## More test runs: run 27 ##
Now we have a transcriptome index file, we can run against this and try different values of -r (as before), but also try a negative value. But now we want to use Bowtie 2 for the mappings as we won't have any spliced reads.

```bash
time bowtie2 -x /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --very-sensitive --no-mixed --no-discordant -1 seqs_100K_1.fastq -2 seqs_100K_2.fastq -p 2 -S run27.sam
```

Because this is Bowtie format, the MAPQ scores will be slightly different (max = 42, not 50), so need to slightly adjust how we extract inner distances:

```bash
cat run27.sam | grep -vE "^@" | awk '$5 == 42 && ($2 == 99 || $2 == 83) {print}'  | grep 101M | cut -f 9 | grep '-' | awk '{print $1,"\t",-$1-202}' > run27_mapping_distances.tsv
```

When plotted in R, this ends up giving a mean inner distance of 34 bp, with a standard deviation of 97 bp. The max inner distance was 298 bp, but this is to be expected because bowtie 2 uses a max fragment length (including reads) of 500. Can change this with -X option, so let's see what happens if you increase -X to 1,000 or 5,000:

```bash
time bowtie2 -x /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --very-sensitive -X 1000 --no-mixed --no-discordant -1 seqs_100K_1.fastq -2 seqs_100K_2.fastq -p 2 -S run28.sam
```

```bash
time bowtie2 -x /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --very-sensitive -X 5000 --no-mixed --no-discordant -1 seqs_100K_1.fastq -2 seqs_100K_2.fastq -p 2 -S run29.sam
```

This latter option only added a couple of extra reads that only had distances marginally greater than 1,000 bp, so -X option doesn't need to go up much above 1,000.

Without changing -X you end up with a mean and standard deviation of 34 and 97 bp. Increasing -X to 1000 seems to stop truncating the data and gives us a mean and standard deviation of 47 and 116 bp. So now let's see the difference of running TopHat again but using `--mate-std-dev` option set to 125 bp. Will also want to see the effect of using the TopHat `--b2-very-sensitive` option which passes the `--very-sensitive` parameter to Bowtie 2. Rather than investigate the average inner difference, just want to compare overall performance (how much output we get):

```
# run30: all defaults (but with Transcriptome settings)
time tophat --no-convert-bam -o run30 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome  -T -M /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq

# run31: default standard deviation (we've probably done this run before)
time tophat --no-convert-bam -o run31 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --no-mixed --no-discordant -T -M /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq

# run32: increase standard deviation to 125
time tophat --no-convert-bam -o run32 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --no-mixed --no-discordant -T -M --mate-std-dev 125 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq

# run33: standard deviation = 125, and use --b2-very-sensitive option
time tophat --no-convert-bam -o run33 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --b2-very-sensitive --no-mixed --no-discordant -T -M --mate-std-dev 125 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq


# run34: standard deviation = 125, -r = 100 and use --b2-very-sensitive option
time tophat --no-convert-bam -o run34 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --b2-very-sensitive --no-mixed --no-discordant -T -M -r 100 --mate-std-dev 125 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq

# run35: standard deviation = 125, -r = 200 and use --b2-very-sensitive option
time tophat --no-convert-bam -o run35 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --b2-very-sensitive --no-mixed --no-discordant -T -M -r 200 --mate-std-dev 125 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq

# run36: standard deviation = 125, -r = 400 and use --b2-very-sensitive option
time tophat --no-convert-bam -o run36 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --b2-very-sensitive --no-mixed --no-discordant -T -M -r 400 --mate-std-dev 125 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq


# run37: standard deviation = 125, -r = 800 and use --b2-very-sensitive option
time tophat --no-convert-bam -o run37 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --b2-very-sensitive --no-mixed --no-discordant -T -M -r 800 --mate-std-dev 125 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq

# run38: standard deviation = 500, -r = 800 and use --b2-very-sensitive option
time tophat --no-convert-bam -o run38 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --b2-very-sensitive --no-mixed --no-discordant -T -M -r 800 --mate-std-dev 500 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq

# run39: standard deviation = 1000, -r = 800 and use --b2-very-sensitive option
time tophat --no-convert-bam -o run39 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --b2-very-sensitive --no-mixed --no-discordant -T -M -r 800 --mate-std-dev 1000 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq

# run40: standard deviation = 1000, -r = 800 no use of --b2-very-sensitive 
time tophat --no-convert-bam -o run40 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --no-mixed --no-discordant -T -M -r 800 --mate-std-dev 1000 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index seqs_100K_1.fastq seqs_100K_2.fastq

```

## Summary of the last set of test runs ##

Best = MAPQ score of 50 and only proper concordant mapped pairs (SAM bitwise flag is 99, 147, 163, or 83)

| Run | Notes                      | Mapping% |  N       | Best     | %Best  |
|:---:|:--------------------------:|:--------:|:--------:|:--------:|:------:|
| 30  | Defaults                   |   73.4%  | 150,235  |  84,996  |  56.6% |
| 31  | No mixed/discordant        |   69.0%  | 139,398  |  83,662  |  60.0% |
| 32  | Std=125                    |   69.0%  | 139,002  | 102,098  |  73.5% |
| 33  | Std=125, sensitive         |   69.0%  | 139,024  | 102,112  |  73.4% |
| 34  | Std=125, r=100, sensitive  |   69.0%  | 138,970  | 105,748  |  76.1% |
| 35  | Std=125, r=200, sensitive  |   69.0%  | 138,848  | 107,112  |  77.1% |
| 36  | Std=125, r=400, sensitive  |   69.0%  | 138,680  | 106,500  |  76.8% |
| 37  | Std=125, r=800, sensitive  |   69.0%  | 139,014  | 111,012  |  79.8% |
| 38  | Std=500, r=800, sensitive  |   69.0%  | 138,792  | 115,770  |  83.4% |
| 39  | Std=1000, r=800, sensitive |   69.0%  | 139,516  | 124,038  |  88.9% |
| 40  | Std=1000, r=800,           |   69.0%  | 139,494  | 124,016  |  88.9% |

On the basis of this, I will proceed using the settings of run 39. This gives the most high quality mapped pairs, even if the value of r may not be 100% appropriate.


# Main run on TopHat

Use `wrapper3.pl` script in `/share/tamu/Analysis/RNA-Seq_FASTQ_files` directory to submit `run_tophat.sh` script to job scheduler. This should create new output directories in `/share/tamu/Analysis/TopHat_output` (one directory for each of the 31 samples). 

Reminder, this run is not going to use the singles files. We may wish to return to these later when looking to improve gene annotations.
