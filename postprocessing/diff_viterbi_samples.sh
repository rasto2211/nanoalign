#!/bin/bash

# This script takes file from HMM sampling. Aligns every sample to Viterbi 
# sequence and processes output from that by get_bwa_stats.py. For some samples
# we might find multiple alignments therefore we append number of the sample in
# front of wvery row. So the first column of output table contains sample
# number. 

file=$1

# First line of input file contains sequence from Viterbi.
# Create separate fasta file for that.
echo ">viterbi" > ${file}_viterbi.fa
head -1 $file >> ${file}_viterbi.fa

bwa index ${file}_viterbi.fa 2>/dev/null

# 2nd line of input is empty. 3rd line contains first sample.
echo "identity" > ${file}.aligned
for line in `sed -n '3,$p' < $file`;
do
  # Create fasta file.
  echo ">sample" > ${file}.fa
  echo $line >> ${file}.fa

  bwa mem ${file}_viterbi.fa ${file}.fa > ${file}.sam 2>/dev/null

  python3 get_bwa_stats.py ${file}.sam |
  python3 choose_greatest_identity_alignment.py >> ${file}.aligned
done

cat ${file}.aligned | ./column_stats.r
