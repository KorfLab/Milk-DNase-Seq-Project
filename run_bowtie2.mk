# Use special .PHONY target for 'clean' and 'all' rules. This is to ensure that 
# 'make clean' will work in the (unlikely?) event that there is a file called clean in the 
# directory without using this special make target, running 'make clean' would fail to run 
# if there was a file called clean. 
.PHONY: clean all


# variable to list all of the source *.fq.gz files that already exist
SOURCES = $(wildcard *_processed.fastq)
MM9_INDEX  = /share/tamu/Data/Bowtie2_indexes/mm9_index
MM10_INDEX = /share/tamu/Data/Bowtie2_indexes/mm10_index

# want to use list of source files to generate target files
# use make's patsubst function to change end of file names
# have two sets for both mm9 and mm10 runs
TARGETS_MM9  = $(patsubst %_processed.fastq, %_mm9.sam,  $(SOURCES))
TARGETS_MM10 = $(patsubst %_processed.fastq, %_mm10.sam, $(SOURCES))



# The main target (all) requires *_processed.fq files
all: $(TARGETS_MM9) $(TARGETS_MM10)

# use make wildcards (%) in rule so that we iterate over all files
# in $(TARGETS_MM9) and $(TARGETS_MM10)
%_mm9.sam: %_processed.fastq
	# $@ refers to the target of the rule (in this case test1_processed.fq)
	# $^ refers to all prequisites of the rule (in this case just test1_scythe.fq)
	# Can also use $< to reflect just the first (or only) pre-requisite
	@echo "Converting $< to $@"
	bowtie2 -x $(MM9_INDEX) -U $< -p 2 --very-sensitive --no-unal --no-hd --omit-sec-seq -S $@

%_mm10.sam: %_processed.fastq
	# $@ refers to the target of the rule (in this case test1_processed.fq)
	# $^ refers to all prequisites of the rule (in this case just test1_scythe.fq)
	# Can also use $< to reflect just the first (or only) pre-requisite
	@echo "Converting $< to $@"
	bowtie2 -x $(MM10_INDEX) -U $< -p 2 --very-sensitive --no-unal --no-hd --omit-sec-seq -S $@


clean:
	@echo "Removing SAM files"
	rm *.sam
