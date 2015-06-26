# Repeatmasking cow genome to find LINE elements

The [RepeatMasker website](http://repeatmasker.org) contains [pre-generated datasets for the bosTau7 version](http://repeatmasker.org/species/bosTau.html) but we want to use bosTau6 (UMD 3.1). So need to repeatmask ourselves.

Running this on saetta (which has the latest version of RepeatMasker installed). We are ultimately only interested in finding LINE elements, so can use some options to exclude classes we are not interested in

## Data prep

```bash
mkdir /Volumes/scratch/keith/Repeatmasking_cow_genome
cp /Korflab/Genomes/Bos_taurus/Bos_taurus_UMD_3.1/bos_taurus.fa.gz .
gunzip bos_taurus.fa.gz
```


## First test using -qq (rush mode)

Use `-nolow`, `-noint`, and `-norna` to turn off masking of some classes of repeats we are not interested in. 

```bash
time RepeatMasker -qq -pa 8 -nolow -noint -norna -gff -species cow bos_taurus.fa
```

This took 35 hours. I renamed the output files to include `_qq`. So now try…


## 2nd test using -s (slow mode)

Should be more sensitive.

```bash
time RepeatMasker -s -pa 8 -nolow -noint -norna -gff -species cow bos_taurus.fa
```

This seemed to crash. Maybe a problem of the -no{*} options? I reported details to RepeatMasker website.


## 3rd test (default mode)

Now trying the regular version (no -q, -qq, or -s):

```bash
time RepeatMasker -pa 8 -gff -species cow bos_taurus.fa
```

Took ~45 hours.


## 4rd test -s (slow mode), again

Back to slow mode, but without the -no{*} options:


```bash
time RepeatMasker -s -pa 8 -gff -species cow bos_taurus.fa
```


## What was the difference between -qq, -s and default mode?

Moving from -qq to default to -s mode gives you more repeats. The total proportion of repeats falls from 47.92% (slow mode) to 46.35% (default mode) to 42.60% (rush mode).




## Also generated some repeats from Repeatmasker.org ##

These are based on bosTau7, but I asked repeatmasker website to find all LINE elements on chr1 (between coordinates 1–2,000,000 bp), with a score > 2,000. This generated 266 sequences which I copied to tamu as a FASTA file.


# On tamu

Make test set and build index:

```bash
mkdir /share/tamu/Analysis/Test/LINE_test
cd /share/tamu/Analysis/Test/LINE_test

head -n 4000000 ../../TopHat_output/L468_C4K66ACXX-7-ID09_1_sequence_processed.fastq > test_1M_1.fastq

head -n 4000000 ../../TopHat_output/L468_C4K66ACXX-7-ID09_2_sequence_processed.fastq > test_1M_2.fastq

bowtie2-build bosTau7_LINE_example_seqs.fa bosTau7_LINEs_index
```


## run 1

Using parameters we've used before (might not be most suitable for this). But supress unaligned reads in output and SAM header:

```
time bowtie2 -x bosTau7_LINEs_index --very-sensitive --no-mixed --no-discordant -1 test_1M_1.fastq -2 test_1M_2.fastq -p 2 --no-unal --no-hd -S run1.sam

1000000 reads; of these:
  1000000 (100.00%) were paired; of these:
    998200 (99.82%) aligned concordantly 0 times
    125 (0.01%) aligned concordantly exactly 1 time
    1675 (0.17%) aligned concordantly >1 times
0.18% overall alignment rate
```


## run 2

Drop requirement for concordancy and mixed. Suppress unaligned reads in output and SAM header:

```
time bowtie2 -x bosTau7_LINEs_index --very-sensitive  -1 test_1M_1.fastq -2 test_1M_2.fastq -p 2 --no-unal --no-hd -S run2.sam

1000000 reads; of these:
  1000000 (100.00%) were paired; of these:
    998200 (99.82%) aligned concordantly 0 times
    125 (0.01%) aligned concordantly exactly 1 time
    1675 (0.17%) aligned concordantly >1 times
    ----
    998200 pairs aligned concordantly 0 times; of these:
      14 (0.00%) aligned discordantly 1 time
    ----
    998186 pairs aligned 0 times concordantly or discordantly; of these:
      1996372 mates make up the pairs; of these:
        1993546 (99.86%) aligned 0 times
        431 (0.02%) aligned exactly 1 time
        2395 (0.12%) aligned >1 times
0.32% overall alignment rate
```


## run 3

Treat as single (unpaired) files:

```bash
time bowtie2 -x bosTau7_LINEs_index --very-sensitive  -U test_1M_1.fastq,test_1M_2.fastq -p 2 --no-unal --no-hd -S run3.sam

2000000 reads; of these:
  2000000 (100.00%) were unpaired; of these:
    1993487 (99.67%) aligned 0 times
    506 (0.03%) aligned exactly 1 time
    6007 (0.30%) aligned >1 times
0.33% overall alignment rate
```


## run 4 and 5

Try against full set of virgin and lactation samples:


```bash
time bowtie2 -x bosTau7_LINEs_index --very-sensitive  -U V468A_C5B2KACXX-5-ID12_1_sequence_processed.fastq,V468A_C5B2KACXX-5-ID12_2_sequence_processed.fastq -p 2 --no-unal --no-hd -S run4.sam

41527888 reads; of these:
  41527888 (100.00%) were unpaired; of these:
    41254361 (99.34%) aligned 0 times
    19313 (0.05%) aligned exactly 1 time
    254214 (0.61%) aligned >1 times
0.66% overall alignment rate

time bowtie2 -x bosTau7_LINEs_index --very-sensitive  -U L468_C4K66ACXX-7-ID09_1_sequence_processed.fastq,L468_C4K66ACXX-7-ID09_2_sequence_processed.fastq -p 2 --no-unal --no-hd -S run5.sam

50050788 reads; of these:
  50050788 (100.00%) were unpaired; of these:
    49891140 (99.68%) aligned 0 times
    12119 (0.02%) aligned exactly 1 time
    147529 (0.29%) aligned >1 times
0.32% overall alignment rate
```

## run 6 and 7

Using more LINEs in database:

```bash
time bowtie2 -x bosTau7_LINEs_index --very-sensitive  -U V468A_C5B2KACXX-5-ID12_1_sequence_processed.fastq,V468A_C5B2KACXX-5-ID12_2_sequence_processed.fastq -p 2 --no-unal --no-hd -S run6.sam

41527888 reads; of these:
  41527888 (100.00%) were unpaired; of these:
    41218588 (99.26%) aligned 0 times
    17340 (0.04%) aligned exactly 1 time
    291960 (0.70%) aligned >1 times
0.74% overall alignment rate

time bowtie2 -x bosTau7_LINEs_index --very-sensitive  -U L468_C4K66ACXX-7-ID09_1_sequence_processed.fastq,L468_C4K66ACXX-7-ID09_2_sequence_processed.fastq -p 2 --no-unal --no-hd -S run7.sam

50050788 reads; of these:
  50050788 (100.00%) were unpaired; of these:
    49868264 (99.64%) aligned 0 times
    9921 (0.02%) aligned exactly 1 time
    172603 (0.34%) aligned >1 times
0.36% overall alignment rate
```