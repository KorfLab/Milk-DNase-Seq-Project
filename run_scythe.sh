#!/bin/bash
# simple bash wrapper around our make file
PATH=$PATH:/home/keith/bin:/share/tamu/bin:/share/tamu/Code
export PATH
cd /share/tamu/Analysis/RNA-Seq_FASTQ_files

make -j 16 -f /share/tamu/Code/run_scythe.mk all

