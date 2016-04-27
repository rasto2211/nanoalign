#!/bin/bash

# Input file contains: ref. read, empty line, lines with seqs. Every seq. is on
# separate line. These seqs. are compared against the ref. read.

file=$1
file_without_ext=${file%.tmp}

# First line of input file contains ref. seq.
echo ">ref" > ${file_without_ext}_ref_seq.fa
head -1 $file >> ${file_without_ext}_ref_seq.fa

bwa index ${file_without_ext}_ref_seq.fa 2>/dev/null

# 2nd line of input is empty. 3rd line contains first seq. which we want to
# compare.
echo "bwa_identity" > ${file_without_ext}_bwa_identity.csv
for line in `sed -n '3,$p' < $file`;
do
  # Create fasta file.
  echo ">seq" > ${file_without_ext}_temp.fa
  echo $line >> ${file_without_ext}_temp.fa

  bwa mem -x ont2d ${file_without_ext}_ref_seq.fa ${file_without_ext}_temp.fa\
  > ${file_without_ext}_temp.sam 2>/dev/null

  # Print the greatest identity or NA.
  python3 get_bwa_stats.py ${file_without_ext}_temp.sam |
  python3 choose_greatest_identity_alignment.py >>\
  ${file_without_ext}_bwa_identity.csv
done

# Delete temp. files.
rm ${file_without_ext}_temp.fa
rm ${file_without_ext}_temp.sam
rm ${file_without_ext}_ref_seq.fa
