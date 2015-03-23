#!/bin/bash
# simple bash script to run TopHat
PATH=$PATH:/home/keith/bin:/share/tamu/bin:/share/tamu/Code
export PATH
cd /share/tamu/Analysis/RNA-Seq_FASTQ_files

# run tophat with two input files ($1 and $2) and put data in one output directory ($3)
echo "Running tophat --no-convert-bam -o $3 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --b2-very-sensitive --no-mixed --no-discordant -T -M -r 800 --mate-std-dev 1000 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index $1 $2"

tophat --no-convert-bam -o $3 --transcriptome-index /share/tamu/Data/Bowtie2_indexes/Ensembl_78_transcriptome --b2-very-sensitive --no-mixed --no-discordant -T -M -r 800 --mate-std-dev 1000 /share/tamu/Data/Bowtie2_indexes/Ensembl-78-index $1 $2