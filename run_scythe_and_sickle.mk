# Use special .PHONY target for 'clean' and 'all' rules. This is to ensure that 
# 'make clean' will work in the (unlikely?) event that there is a file called clean in the 
# directory without using this special make target, running 'make clean' would fail to run 
# if there was a file called clean. 
.PHONY: clean all

# need variable (macro)( name for adapter sequences
ADAPTER_SEQUENCES = /share/tamu/Data/DNase-Seq/Mouse/adapter_sequences.fasta

# variable to list all of the source *.fq.gz files that already exist
SOURCES = $(wildcard *.fastq.gz)

# want to use list of source files to generate target files
# use make's patsubst function to change end of file names
TARGETS = $(patsubst %.fastq.gz, %_processed.fastq, $(SOURCES))



# The main target (all) requires *_processed.fq files
all: $(TARGETS)


# use make wildcards (%) in rule so that we iterate over all files
# in $(TARGETS)
%_processed.fastq: %_scythe.fastq
	# $@ refers to the target of the rule (in this case test1_processed.fq)
	# $^ refers to all prequisites of the rule (in this case just test1_scythe.fq)
	# Can also use $< to reflect just the first (or only) pre-requisite
	# so following two lines are equivalent
	# sickle se -f test1_scythe.fq -o test1_processed.fq -t sanger
	sickle se -f $< -o $@ -t sanger
	# once we have run sickle, can remove the scythe file
	@echo "Removing intermediate Scythe file"
	rm $<

# this will be the first main rule that does anything
# run Scythe on each FASTQ file to make a scythe output file
%_scythe.fastq: %.fastq.gz
	scythe -a $(ADAPTER_SEQUENCES) -o $@ $<


# If there were new adapter sequences, we would need to update Scythe,
# can make the adapters part of a prerequisite for the same rule that is
# used to run Scythe. When make has multiple rules with the same name
# it takes the union of them. This allows us to achieve what we want
#test1_scythe.fq: $(ADAPTER_SEQUENCES)
#test2_scythe.fq: $(ADAPTER_SEQUENCES)

# The above works if we no how many files there are
# but if you have unknown number you'd ideally like to use a make wildcard
# character in this rule. But if make sees two rules using patterns,
# it will only execute the first one. So need to make false dependencies ruleÉ
# We know that FQ files will exist so we can make them the target for checking 
# that adapter_sequences.fasta file exists
*.fastq : $(ADAPTER_SEQUENCES)
	touch $@

clean:
	@echo "Removing processed FASTQ files"
	rm *_processed.fq
