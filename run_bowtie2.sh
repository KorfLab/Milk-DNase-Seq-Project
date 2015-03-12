#!/bin/bash
# simple bash wrapper around our make file
PATH=$PATH:/home/keith/bin:/share/tamu/bin:/share/tamu/Code
export PATH
cd /share/tamu/Analysis/All_FASTQ_files

make -j 4 -f /share/tamu/Code/run_bowtie2.mk all

