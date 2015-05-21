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
qsub -S /bin/bash -pe threaded 1 -l h_vmem=8G -M keith\@bradnam.co -m be -N krb_tophat_run /share/tamu/Code/run_tophat.sh L478_C4K66ACXX-7-ID10_1_sequence_processed.fastq L478_C4K66ACXX-7-ID10_2_sequence_processed.srr /share/tamu/Analysis/TopHat_output/L478_C4K66ACXX-7-ID10
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
cat test_1M.sam | w

samtools  view -q 50 -f 0x2 test_1M.sam


cat test_1M.sam | awk '$5 == 50 && ($2 == 163 || $2 == 147 || $2 == 83 || $2 == 99) {print}' > test_1M_HQ.sam

cat test_5M.sam | awk '$5 == 50 && ($2 == 163 || $2 == 147 || $2 == 83 || $2 == 99) {print}' > test_5M_HQ.sam

cat test_25M.sam | awk '$5 == 50 && ($2 == 163 || $2 == 147 || $2 == 83 || $2 == 99) {print}' > test_25M_HQ.sam

```

# Testing HTSeq

Need to use `-r pos` to tell htseq-count that our data is already sorted by position. The default GFF feature to map to is exon (`-t` option) which should suffice as should the default value of 'gene_id' for the `-i` option (specifying the attribute field). First runs will test difference of using `-s no` to specify that this is not a stranded assay (default is yes). Using a full SAM file (>31 million reads):

```bash

htseq-count -r pos -m intersection-nonempty L468_C4K66ACXX-7-ID09/accepted_hits.sam ../Bowtie2_indexes/Ensembl_78_transcriptome.gff > test_htseq_run1.txt

htseq-count -r pos -m intersection-nonempty -s no L468_C4K66ACXX-7-ID09/accepted_hits.sam ../Bowtie2_indexes/Ensembl_78_transcriptome.gff > test_htseq_run2.txt

tail -n 5 test_htseq_run*
==> test_htseq_run1.txt <==
__no_feature    15848564
__ambiguous     10
__too_low_aQual 0
__not_aligned   0
__alignment_not_unique  0

==> test_htseq_run2.txt <==
__no_feature    10868
__ambiguous     2769
__too_low_aQual 0
__not_aligned   0
__alignment_not_unique  0

```


We definitely need to use `-s no` otherwise it doesn't count reads that are on the opposite strand to an exon feature.

# Running HTseq

Wrote a simple wrapper Perl script to first filter SAM files (using awk approach shown above) and then run htseq-count. Final output files are in `/share/tamu/Analysis/Filter_SAM_files`

I then wrote another Perl wrapper script `/share/tamu/Code/join_htseq_count_files.pl` that uses the Unix `join` command to join all of the counts files into one final file:

```bash
join_htseq_count_files.pl "*-ID*.txt"
mv tmp_join_file.txt count_data_for_all_samples.tsv
```

This file includes error counts in last few rows which we want to remove:

```bash
grep -v __ count_data_for_all_samples.tsv  > tmp; mv tmp count_data_for_all_samples.tsv
```

Could also filter each file separately:

```bash
for f in *.txt; do grep -v __ $f > tmp; mv -f tmp $f; done
```

# DESseq analysis

Installed latest DESeq2 package in R using following R commands:

```R
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
```

Following DESeq2 information in these online guides and manuals:

1. <http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/>
2. <http://bioconductor.fmrp.usp.br/packages/2.14/bioc/vignettes/DESeq2/inst/doc/beginner.pdf>
3. <http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf>
4. <http://www.sthda.com/english/wiki/rna-seq-differential-expression-work-flow-using-deseq2>
5. <http://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf>

Borrowing from other people's examples, it seems that you can combine multiple htseq count files within the DESeq2 package (obviating the need for some of my earlier steps). So I copied count data from individual files into a subdirectory (`HTSeq_count_files/`). The R script `DESeq2_analysis.R` was used to generate log-odds plot, gene specific expression plots, and the final CSV output files.

For an overview of what DESeq revealed, see the `summary_for_monique.md` file.



# DE analysis using limma/voom

Following limma information in this user guide <http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf> and writing code in [limma_analysis](limma_analysis.R) script.



# Applying Dilution Effect

I ended up rewriting and simplifying a lot of Kristen's code. The forked repository on GitHub is at <https://github.com/kbradnam/DilutionEffect>.

Having split my main counts file into virgin and lactation samples, I ran her code like so:

```bash
./dilution.pl virgin_counts.tsv lactation_counts.tsv

# Process counts file to determine mean expression values for each replicate

# Testing for difference in mean expression levels and goodness of fit between both conditions
  Running Mann-Whitney-Wilcoxon test to compare means: P =  0.1201827
  Running Kolmogorov-Smirnov test to compare distributions: P =  0.00000000001608391
  Underlying distributions are significantly different but mean expression is not.
  There may be some value in applying dilution adjustment.

# Applying dilution factor to virgin_counts.tsv
  Calculating thresholds using quantile = 0.9995
  Determining high abundance genes
  Writing the adjusted and pseudocounted data to files

# Applying dilution factor to lactation_counts.tsv
  Calculating thresholds using quantile = 0.9995
  Determining high abundance genes
  Writing the adjusted and pseudocounted data to files
```

So the mean expression in virgin vs lactation is not significantly different, but the overall distributions are. Following Danielle's suggestion, I normalized all of the raw counts (dividing by the total expression in each replicate) by using the --normalize option:

```bash

./dilution.pl --normalize virgin_counts.tsv lactation_counts.tsv

# Testing for difference in mean expression levels and goodness of fit between both conditions
# Using normalized count data
  Running Mann-Whitney-Wilcoxon test to compare means: P =  1.218537e-16
  Running Kolmogorov-Smirnov test to compare distributions: P =  0
  Average expression and underlying distributions are statistically different.
  The dilution adjustment should be applied.
```

This now reports both tests as being significant.

# Assessing effect of DESeq2 on zero counted genes #

Kristen's code pseudocounts genes so that you don't end up with situations with all replicates having counts of zero (leading to a mean of zero which presumably can't be used in DE analysis). This could throw away potentially useful data.

Need to check whether this might be a common situation in our data. I extracted all of the genes with zero counts across all replicates in both the virgin and lactation conditions:

```bash
cat count_data_for_all_samples.tsv  | cut -f 1,18-31 | tr '\t' ',' | grep "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0" | sed 's/,.*//'  > virgin_genes_with_zero_counts.txt

cat count_data_for_all_samples.tsv  | cut -f 1,2-17 | tr '\t' ',' | grep "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0" | sed 's/,.*//'  > lactation_genes_with_zero_counts.txt

```

Now to see how many genes there are with zero counts in each set, and what the overlap is between the two (i.e. how many genes never show any expression in any replicate of either condition?).

```bash
wc -l virgin_genes_with_zero_counts.txt lactation_genes_with_zero_counts.txt
     703 virgin_genes_with_zero_counts.txt
     719 lactation_genes_with_zero_counts.txt

grep -cf virgin_genes_with_zero_counts.txt lactation_genes_with_zero_counts.txt
654

grep -cf lactation_genes_with_zero_counts.txt virgin_genes_with_zero_counts.txt
654     
```

So the majority of cases (654/703 or 654/719) are where there are zero counts in all replicates. 

```bash

grep -vf lactation_genes_with_zero_counts.txt virgin_genes_with_zero_counts.txt > only_expressed_in_lactation.txt

grep -vf virgin_genes_with_zero_counts.txt lactation_genes_with_zero_counts.txt > only_expressed_in_virgin.txt
```


We can then see more details about these genes which are expressed in only one condition (output not included here):
```bash

grep -f only_expressed_in_virgin.txt count_data_for_all_samples.tsv | cut -f 1,18-31 
grep -f only_expressed_in_lactation.txt count_data_for_all_samples.tsv | cut -f 1-17 

grep -f only_expressed_in_virgin.txt deseq2_output_virgin_vs_lactation.csv
grep -f only_expressed_in_lactation.txt deseq2_output_virgin_vs_lactation.csv

```

All of this reveals that:

+ There are 703 genes in virgin replicates where all replicates have zero counts
+ There are 719 genes in lactation replicates where all replicates have zero counts
+ There are 654 genes with zero counts in *both* conditions. These end up with NA values in the final DESeq2 analysis and are not of any interest to us.
+ Genes with zero counts in *one* condition are still analyzed by DESeq2
+ Only 1 of 49 genes expressed in lactation (with zero counts in virgin) is classified as upregulated (log2fold change > 2): ENSBTAG00000010730
+ 2 of these 49 genes end up with NA values. Not obvious why. They have counts in only one replicate, but this is true of other genes. Maybe because of high variance (the counts in question are 12 and 32).
+ 65 genes are only expressed in virgin condition, and only 1 (ENSBTAG00000012390) has a log2fold change greater than -2 (i.e. downregulated in lactation)
+ In nearly all cases, the lactation replicates that have non-zero counts  (where virgin counts are all zero) show very low counts (mostly 0-2). Only one replicate stands out having much higher counts for some of the 49 genes (L581), with counts as high as 64 for one gene.


## Notes on normalization

From an [online guide](http://chipster.csc.fi/manual/deseq2.html):

>"DESeq2 performs an internal normalization where geometric mean is calculated for each gene across all samples. The counts for a gene in each sample is then divided by this mean. The median of these ratios in a sample is the size factor for that sample. This procedure corrects for library size and RNA composition bias, which can arise for example when only a small number of genes are very highly expressed in one experiment condition but not in the other. "

From the [original DESeq paper](http://www.biomedcentral.com/content/pdf/gb-2010-11-10-r106.pdf):

>"The total number of reads, Σi kij, may seem to be a good measure of sequencing depth and hence a reasonable choice for sj. Experience with real data, however, shows this not always to be the case, because a few highly and differentially expressed genes may have strong influence on the total read count, causing the ratio of total read counts not to be a good estimate for the ratio of expected counts.
>
>Hence, to estimate the size factors, we take the median of the ratios of observed counts."

From the manual: 

>"The count values must be raw counts of sequencing reads. This is important for DESeq2’s statistical model [1] to hold, as only the actual counts allow assessing the measurement precision correctly. Hence, please do not supply other quantities, such as (rounded) normalized counts, or counts of covered base pairs – this will only lead to nonsensical results."

>"A simple function for making this plot is plotCounts, which normalizes counts by sequencing depth and adds a pseudocount of 1 to allow for log scale plotting"

From a paper [comparing normalization methods](http://bib.oxfordjournals.org/content/14/6/671.full#ref-14)

>"This normalization method [14] is included in the DESeq Bioconductor package (version 1.6.0) [14] and is based on the hypothesis that most genes are not DE. A DESeq scaling factor for a given lane is computed as the median of the ratio, for each gene, of its read count over its geometric mean across all lanes. The underlying idea is that non-DE genes should have similar read counts across samples, leading to a ratio of 1. Assuming most genes are not DE, the median of this ratio for the lane provides an estimate of the correction factor that should be applied to all read counts of this lane to fulfill the hypothesis. By calling the estimateSizeFactors() and sizeFactors() functions in the DESeq Bioconductor package, this factor is computed for each lane, and raw read counts are divided by the factor associated with their sequencing lane."


# Installing cufflinks to calculate RPKM/FPKM values

Will install v.2.2.1 On tamu:

```bash
mkdir /share/tamu/Packages/Cufflinks
cd /share/tamu/Packages/Cufflinks
curl -O http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar -xvzf cufflinks-2.2.1.Linux_x86_64.tar.gz
ln -s cufflinks-2.2.1.Linux_x86_64 current
cd /share/tamu/bin/
ln -s ../Packages/Cufflinks/current/cufflinks
ln -s ../Packages/Cufflinks/current/cuffdiff
ln -s ../Packages/Cufflinks/current/cuffcompare
ln -s ../Packages/Cufflinks/current/cuffmerge
ln -s ../Packages/Cufflinks/current/cuffnorm
ln -s ../Packages/Cufflinks/current/cuffquant
ln -s ../Packages/Cufflinks/current/gffread
ln -s ../Packages/Cufflinks/current/gtf_to_sam
```

# Testing cufflinks to calculate FPKM

First set up some test data (10 million lines of a SAM file):

```bash
mkdir /share/tamu/Analysis/Test/Cufflinks_test
cd /share/tamu/Analysis/Test/Cufflinks_test
head -n 1000000 /share/tamu/Analysis/Filter_SAM_files/L468_C4K66ACXX-7-ID09_hq.sam > L468_10M.sam
```


## Test 1: just see if it works at all

```bash
time cufflinks -o run1 L468_10M.sam
```

Took 1 minute 56 seconds.


## Test 2: does -p 2 speed things up?

```bash
time cufflinks -o run2 -p 2 L468_10M.sam
```

Down to 1 minute 13 seconds.


## Test 3: specifying a GTF file

```bash
ln -s /share/tamu/Data/Genomes/Cow/bosTau6/Ensembl-78-annotations/Bos_taurus.UMD3.1.78.gtf
time cufflinks -o run3 -p 2 --GTF Bos_taurus.UMD3.1.78.gtf L468_10M.sam
```

Up to 3 minutes 20 seconds.


## Test 4: Include --multi-read-correct option

"Tells Cufflinks to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome"

```bash
time cufflinks -o run4 -p 2 --multi-read-correct --GTF Bos_taurus.UMD3.1.78.gtf L468_10M.sam
```

Up to 4 minutes 38 seconds.



## Test 5: Include --max-bundle-frags option

"Sets the maximum number of fragments a locus may have before being skipped" — Danielle had previously increased the default value from 1 million to 500 million.

```bash
time cufflinks -o run5 -p 2 --multi-read-correct --max-bundle-frags 500000000 --GTF Bos_taurus.UMD3.1.78.gtf L468_10M.sam
```

Slightly up at 4 minutes 45 seconds.



## Test 6: Try –-frag-bias-correct option

"Providing Cufflinks with a multifasta file via this option instructs it to run our new bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates." 

```bash
ln -s ../../../Data/Genomes/Cow/bosTau6/Ensembl-78-genome/Bos_taurus.UMD3.1.dna.toplevel.fa
time cufflinks -o run6 -p 2 --multi-read-correct --max-bundle-frags 500000000 --GTF Bos_taurus.UMD3.1.78.gtf --frag-bias-correct Bos_taurus.UMD3.1.dna.toplevel.fa L468_10M.sam
```

Up to 6 minutes 10 seconds.


I used these settings on Isner with a [patched version of Cufflinks](https://groups.google.com/forum/#!topic/tuxedo-tools-users/UzLCJhj3lUE) because the official version seemed to keep on hanging at a certain point in each input file. The patched version wouldn't run on merlot.

I copied the results directories to `/share/tamu/Analysis/RNA-Seq_Cufflinks_output`


# Mapping RNA-seq reads to entire genome

This will enable the following:

1. Assess existing annotation of bosTau6 genome
2. Look for novel genes and gene variants missed by existing annotation
3. Look to see whether LINE elements are transcribed


Already have bowtie2 index for genome from earlier step (`/share/tamu/Data/Bowtie2_indexes/bosTau6_index.*`)

## Test run

Trying `--time` to 'print the wall-clock time required to load the index files and align the reads'. Other options are used as before, generally being very conservative by running in `--very-sensitive` mode with `--no-mixed` and `--no-discordant` options which will end up throwing away data where only one read from a pair matches well.

Will start with one pair of FASTQ files to get a feel for how long this will take.

```bash
time bowtie2 -x /share/tamu/Data/Bowtie2_indexes/bosTau6 --very-sensitive --no-mixed --no-discordant -X 1000 --time -1 L468_C4K66ACXX-7-ID09_1_sequence_processed.fastq -2 L468_C4K66ACXX-7-ID09_2_sequence_processed.fastq -S L468.sam
```











