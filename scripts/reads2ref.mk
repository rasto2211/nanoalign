# Usage: make -j8 -f reads2ref.mk INPUT=~/2d_reads_move2

FAST5 = $(wildcard $(INPUT)/*.fast5)
ALIGNED_READS = $(FAST5:.fast5=_aligned.fa)
REF_SEQ = ~/ecoliref.fas

all: $(REF_SEQ).bwt $(ALIGNED_READS) $(INPUT)/identity.csv

$(REF_SEQ).bwt: $(REF_SEQ)
	bwa index $(REF_SEQ)

%.fasta: %.fast5
	# Take template from the fasta file.
	poretools fasta $*.fast5 | tail -2 > $*.fasta

%.sam: %.fasta
	~/bwa/bwa mem -x ont2d $(REF_SEQ) $*.fasta > $*.sam 2>/dev/null

%_aligned.fa: %.sam
	python3 correct_read_from_sam.py $*.sam > $*_aligned.fa

$(INPUT)/identity.csv: $(ALIGNED_READS)
	cat $(ALIGNED_READS) | grep '^>.*' | cut -d' ' -f3 > $(INPUT)/identity.csv
	touch $(INPUT)/nieco.csv
