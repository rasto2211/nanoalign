FAST5 = $(wildcard $(INPUT)/*.fast5)
ALIGNED_READS = $(FAST5:.fast5=_aligned.fa)
REF_SEQ = ecoliref.fas

all: $(REF_SEQ).bwt $(ALIGNED_READS)

$(REF_SEQ).bwt: $(REF_SEQ)
	source ~/env/bin/activate
	bwa index $(REF_SEQ)

%.sam: %.fast5
	~/bwa/bwa mem -x ont2d $(REF_SEQ) $* > $*.sam 2>/dev/null

%_aligned.fa: %.sam
	python3 correct_read_from_sam.py $*.sam > $*_aligned.fa
