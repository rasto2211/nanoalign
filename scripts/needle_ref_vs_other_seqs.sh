#!/bin/bash

# Input file contains: ref. read, empty line, lines with seqs. Every seq. is on
# separate line. These seqs. are compared against the ref. read.

file=$1
file_without_ext=${file%.needle_tmp}

# First line of input file contains ref. seq.
echo ">ref" > ${file_without_ext}_ref_seq_needle.fa
head -1 $file >> ${file_without_ext}_ref_seq_needle.fa

# 2nd line of input is empty. 3rd line contains first seq. which we want to
# compare.
echo "needle_identity" > ${file_without_ext}_needle_identity.csv
for line in `sed -n '3,$p' < $file`;
do
  # Create fasta file.
  echo ">seq" > ${file_without_ext}_temp_needle.fa
  echo $line >> ${file_without_ext}_temp_needle.fa

  ./run_needle.sh ${file_without_ext}_ref_seq_needle.fa \
  ${file_without_ext}_temp_needle.fa >> ${file_without_ext}_needle_identity.csv
done

# Delete temp. files.
rm ${file_without_ext}_temp_needle.fa
rm ${file_without_ext}_ref_seq_needle.fa
