#!/bin/bash
# simple bash wrapper around our make file
PATH=$PATH:/home/keith/bin:/share/tamu/bin:/share/tamu/Code
export PATH
cd /share/tamu/Analysis/All_FASTQ_files

make -j  -f /share/tamu/Code/run_scythe_and_sickle.mk all

