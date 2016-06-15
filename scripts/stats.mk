# This script runs the whole pipeline for analizing Viterbi seq. and samples.
# Usage: make -j8 -f stats.mk INPUT=~/2d_reads_move2

REF_SEQ = ~/ecoliref.fas 

FAST5 = $(wildcard $(INPUT)/*.fast5)
SAMPLES = $(wildcard $(INPUT)/*.samples)
ALIGNED_READS = $(FAST5:.fast5=_aligned.fa)
# Sample with part of ref. seq. corresponding to this read.
SAMPLES_REF = $(SAMPLES:.samples=.samples_ref) 
INTERSECTION_CSV = $(SAMPLES_REF:.samples_ref=_intersection.csv)
ALL_IN_ONE_PLOTS = $(INTERSECTION_CSV:_intersection.csv=_all_in_one.pdf)
SEPARATE_PLOTS = $(INTERSECTION_CSV:_intersection.csv=_separate.pdf)
BWA_CSV = $(SAMPLES:.samples=_bwa.csv) 
NEEDLE_IDENTITY = $(SAMPLES:.samples=_needle_identity.csv) 
BASELINES_CSV = $(SAMPLES_REF:.samples_ref=_baselines.csv)

.SECONDARY:
# Prevents deletion of intermediate targets in chained rules.

.DELETE_ON_ERROR:
# Delete targets if a rule fails.

all: $(ALIGNED_READS) $(SAMPLES_REF) $(INTERSECTION_CSV)\
	$(ALL_IN_ONE_PLOTS) $(SEPARATE_PLOTS)\
	plots_all_in_one.pdf plots_separate.pdf $(BWA_CSV)\
	compound_intersection.pdf all_bwa.csv bwa_identity.pdf 
	# needle_identity.pdf -- too much time consuming and does not give good 
	# results.

# Index ref.
$(REF_SEQ).bwt: $(REF_SEQ)
	bwa index $(REF_SEQ) 2>/dev/null

%.samples_no_pipes: %.samples
	cat $*.samples | tr -d '|' > $*.samples_no_pipes

# Extract metrichor basecall to .fast5
%.fasta: %.fast5
	# Take template from the fasta file.
	poretools fasta $*.fast5 | tail -2 > $*.fasta

# Construct BWA-MEM alignment.
%.sam: %.fasta $(REF_SEQ).bwt 
	~/bwa/bwa mem -x ont2d $(REF_SEQ) $*.fasta > $*.sam 2>/dev/null

# Extract best alignment from SAM. Output contains name seq.
# and part of ref. seq. that the read was aligned to.
%_aligned.fa: %.sam
	python3 correct_read_from_sam.py $*.sam > $*_aligned.fa

# Output contains part of ref. seq. aligned to the read, empty line and samples.
%.samples_ref: %.samples_no_pipes %_aligned.fa
	# Copy ref. read.
	tail -1 $*_aligned.fa > $*.samples_ref
	sed -n '2,$$p' < $*.samples_no_pipes >> $*.samples_ref

# Compute size of intersection of kmers between ref. read and samples.
%_intersection.csv: %.samples_ref
	../src/kmers_intersection_samples_main --samples_file=$*.samples_ref \
	--k_low=9 --k_upper=30 --logtostderr > $*_intersection.csv

# Collect sequences for Viterbi and Metrichor for a particular read.
%_baselines.seqs: %.samples_no_pipes %.fasta %_aligned.fa
	# Copy ref. read.
	tail -1 $*_aligned.fa > $*_baselines.seqs
	# Empty line
	echo "" >> $*_baselines.seqs
	# Copy Viterbi sequence.
	head -1 $*.samples_no_pipes >> $*_baselines.seqs 
	# Copy metrichor basecalled sequence.
	tail -1 $*.fasta >> $*_baselines.seqs

# Compare baselines with ref. seq.
%_baselines.csv: %_baselines.seqs
	../src/kmers_intersection_seqs_main --seqs_file=$*_baselines.seqs \
	--k_low=9 --k_upper=30 --logtostderr > $*_baselines.csv

# These graphs contain plots for all k.
%_all_in_one.pdf: %_baselines.csv %_intersection.csv
	Rscript plot_all_in_one.r $*_intersection.csv $*_baselines.csv $*

# Merge ALL_IN_ONE_PLOTS to one pdf.
plots_all_in_one.pdf: $(ALL_IN_ONE_PLOTS)
	pdftk $(ALL_IN_ONE_PLOTS) cat output plots_all_in_one.pdf

# Separate graph for every k.
%_separate.pdf: %_baselines.csv %_intersection.csv
	Rscript plot_separate.r $*_intersection.csv $*_baselines.csv $*

# Merge SEPARATE_PLOTS to one pdf.
plots_separate.pdf: $(SEPARATE_PLOTS)
	pdftk $(SEPARATE_PLOTS) cat output plots_separate.pdf

# Collect Needle identities for Viterbi, Metrichor and samples per read.
%_needle_identity.csv: %.samples_ref %.fasta %.samples_no_pipes
	# Take ref. read and samples.
	cat $*.samples_ref > $*.needle_tmp
	# Take Viterbi seq.
	head -1 $*.samples_no_pipes >> $*.needle_tmp
	# Take Metrichor seq.
	tail -1 $*.fasta >> $*.needle_tmp
	./needle_ref_vs_other_seqs.sh $*.needle_tmp

# Collect BWA stats for Viterbi, Metrichor and samples per read.
%_bwa.csv: %.samples_ref %.fasta %.samples_no_pipes
	# Take ref. read and samples.
	cat $*.samples_ref > $*.tmp
	# Take Viterbi seq.
	head -1 $*.samples_no_pipes >> $*.tmp
	# Take Metrichor seq.
	tail -1 $*.fasta >> $*.tmp
	# Produce $*_bwa.csv
	./bwa_ref_vs_other_seqs.sh $*.tmp |\
	python3 csv_with_read_stats.py 1 > $*_bwa.csv

all_bwa.csv: $(BWA_CSV)
	echo $(BWA_CSV) | tr ' ' '\n' | python3 merge_csv.py > all_bwa.csv

# Produce BWA identity box plot
bwa_identity.pdf: all_bwa.csv
	python3 sample_reads.py 40 all_bwa.csv
	Rscript identity_box_plot.r id_samples.csv bwa_identity_clip.pdf
	Rscript identity_box_plot.r id_clip_samples.csv bwa_identity_no_clip.pdf
	pdftk bwa_identity_clip.pdf bwa_identity_no_clip.pdf \
	cat output bwa_identity.pdf
	rm bwa_identity_clip.pdf
	rm bwa_identity_no_clip.pdf

needle_identity.pdf: $(NEEDLE_IDENTITY)
	echo $(NEEDLE_IDENTITY) | tr ' ' '\n' |\
	python3 sample_reads_identities.py 35 >\
	$(INPUT)/needle_identities_samples.csv
	Rscript identity_box_plot.r $(INPUT)/needle_identities_samples.csv \
	needle_identity.pdf

# Merge all $(INTERSECTION_CSV) files into one CSV.
$(INPUT)/compound_samples_intersection.csv: $(INTERSECTION_CSV)
	echo $(INTERSECTION_CSV) | tr ' ' '\n' |\
	python3 create_compound_csv.py >\
	$(INPUT)/compound_samples_intersection.csv

# Merge all $(BASELINES_CSV) file into one CSV.
$(INPUT)/compound_baselines_intersection.csv: $(BASELINES_CSV)
	echo $(BASELINES_CSV) | tr ' ' '\n' | python3 create_compound_csv.py >\
	$(INPUT)/compound_baselines_intersection.csv

# Create compound graphs.
compound_intersection.pdf: $(INPUT)/compound_samples_intersection.csv\
$(INPUT)/compound_baselines_intersection.csv
	Rscript plot_compound_intersection.r \
	$(INPUT)/compound_samples_intersection.csv \
	$(INPUT)/compound_baselines_intersection.csv
