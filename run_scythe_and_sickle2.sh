#!/bin/bash
# simple bash wrapper around our make file
PATH=$PATH:/home/keith/bin:/share/tamu/bin:/share/tamu/Code
export PATH
cd /share/tamu/Analysis/RNA-Seq_FASTQ_files

# first scythe file
echo "Running scythe -a /share/tamu/Data/DNase-Seq/Mouse/adapter_sequences.fasta -o $3 $1"
scythe -a /share/tamu/Data/DNase-Seq/Mouse/adapter_sequences.fasta -o $3 $1

# second scythe file
echo "Running scythe -a /share/tamu/Data/DNase-Seq/Mouse/adapter_sequences.fasta -o $4 $2"
scythe -a /share/tamu/Data/DNase-Seq/Mouse/adapter_sequences.fasta -o $4 $2


# paired end sickle run
echo "Running sickle pe -f $3 -r $4 -t sanger -o $5 -p $6 -s $7"
sickle pe -f $3 -r $4 -t sanger -o $5 -p $6 -s $7

