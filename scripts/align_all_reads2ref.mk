# Example of running this makefile:
# make -f align_all_reads2ref.mk INPUT=~/ecoli

REF_SEQ = ~/ecoliref.fas 

FAST5 = $(wildcard $(INPUT)/*.fast5)
SAM = $(FAST5:.fast5=.sam)

all: $(SAM)

# Index ref.
$(REF_SEQ).bwt: $(REF_SEQ)
	bwa index $(REF_SEQ) 2>/dev/null

# Extract metrichor basecall to .fast5
%.fasta: %.fast5
	# Take template from the fasta file.
	poretools fasta $*.fast5 | tail -2 > $*.fasta

# Construct BWA-MEM alignment.
%.sam: $(REF_SEQ).bwt %.fasta  
	~/bwa/bwa mem -x ont2d $(REF_SEQ) $*.fasta > $*.sam 2>/dev/null
