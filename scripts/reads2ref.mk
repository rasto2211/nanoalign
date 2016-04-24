# Usage: make -j8 -f reads2ref.mk INPUT=~/2d_reads_move2

REF_SEQ = ~/ecoliref.fas 

FAST5 = $(wildcard $(INPUT)/*.fast5)
SAMPLES = $(wildcard $(INPUT)/*.samples)
ALIGNED_READS = $(FAST5:.fast5=_aligned.fa)
# Sample with part of ref. seq. corresponding to this read.
SAMPLES_REF = $(SAMPLES:.samples=.samples_ref) 
CSV = $(SAMPLES_REF:.samples_ref=_intersection.csv)
ALL_IN_ONE_PLOTS = $(CSV:_intersection.csv=_all_in_one.pdf)
SEPARATE_PLOTS = $(CSV:_intersection.csv=_separate.pdf)

.SECONDARY:
# Prevents deletion of intermediate targets in chained rules.

.DELETE_ON_ERROR:
# Delete targets if a rule fails.

all: $(ALIGNED_READS) $(SAMPLES_REF) $(INPUT)/identity.csv $(CSV) $(ALL_IN_ONE_PLOTS) $(SEPARATE_PLOTS) plots_all_in_one.pdf plots_separate.pdf

$(REF_SEQ).bwt: $(REF_SEQ)
	bwa index $(REF_SEQ) 2>/dev/null

%.fasta: %.fast5
	# Take template from the fasta file.
	poretools fasta $*.fast5 | tail -2 > $*.fasta

%.sam: %.fasta $(REF_SEQ).bwt 
	~/bwa/bwa mem -x ont2d $(REF_SEQ) $*.fasta > $*.sam 2>/dev/null

%_aligned.fa: %.sam
	python3 correct_read_from_sam.py $*.sam > $*_aligned.fa

%.samples_ref: %.samples %_aligned.fa
	# Copy ref. read.
	tail -1 $*_aligned.fa > $*.samples_ref
	sed -n '2,$$p' < $*.samples >> $*.samples_ref

%_intersection.csv: %.samples_ref
	../src/kmers_intersection_samples_main --samples_file=$*.samples_ref \
	--k_low=9 --k_upper=30 --step=15 > $*_intersection.csv

%_baselines.seqs: %.samples %.fasta %_aligned.fa
	# Copy ref. read.
	tail -1 $*_aligned.fa > $*_baselines.seqs
	# Empty line
	echo "" >> $*_baselines.seqs
	# Copy Viterbi sequence.
	head -1 $*.samples >> $*_baselines.seqs 
	# Copy metrichor basecalled sequence.
	tail -1 $*.fasta >> $*_baselines.seqs

# Compare baselines with ref. seq.
%_baselines.csv: %_baselines.seqs
	../src/kmers_intersection_seqs_main --seqs_file=$*_baselines.seqs \
	--k_low=9 --k_upper=30 > $*_baselines.csv

%_all_in_one.pdf: %_baselines.csv %_intersection.csv
	Rscript plot_all_in_one.r $*_intersection.csv $*_baselines.csv $*

plots_all_in_one.pdf: $(ALL_IN_ONE_PLOTS)
	pdftk $(ALL_IN_ONE_PLOTS) cat output plots_all_in_one.pdf

%_separate.pdf: %_baselines.csv %_intersection.csv
	Rscript plot_separate.r $*_intersection.csv $*_baselines.csv $*

plots_separate.pdf: $(SEPARATE_PLOTS)
	pdftk $(SEPARATE_PLOTS) cat output plots_separate.pdf

# Create .csv with all identity percentages after running bwa-mem.
$(INPUT)/identity.csv: $(ALIGNED_READS)
	grep --no-filename '^>.*' $(ALIGNED_READS) | sed -r 's/^>(.*)/\1/' > $(INPUT)/identity.csv
