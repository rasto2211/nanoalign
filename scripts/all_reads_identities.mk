# Example of runnint this makefile:
# make -j4 -f all_reads_identities.mk INPUT:=~/samples

SAMPLES = $(wildcard $(INPUT)/*.samples)

all: statistics.tsv

statistics.tsv: $(SAMPLES:.samples=.tsv)
	echo 'sd\tmean\tmedian\tNAs\t25%\t75%\tmin\tmax' > statistics.tsv
	cat $^ >> statistics.tsv
	rm $^

%.tsv : %.samples 
	./bwa_ref_vs_other_seqs.sh $*.samples
	cat $*_bwa_identity.csv | ./column_stats.r > $*.tsv
