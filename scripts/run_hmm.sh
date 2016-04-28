#!/bin/bash

# Takes file with list of paths of reads and path of working folder in which all
# files produced by HMM will be placed.

training_set_percent=0.7
move_threshold=2
pseudocount=1
samples=250

input_reads_list=$1
working_folder=$2
scripts_folder=~/nanopore-read-align/scripts
strand=template
template_strand=true

mkdir -p $working_folder
cd $working_folder;
mkdir -p hmm;
cd hmm;

while read line;
do
    ln -s $line ../
done < ${input_reads_list};

echo 'Finished creating symlinks. Splitting data...';

python3 ${scripts_folder}/split_data.py \
--file_with_list=${input_reads_list} \
--strand=$strand \
--training_set_percent=$training_set_percent

echo 'Data set was successfully split. Trainning HMM...';

${scripts_folder}/../src/train_move_hmm_main \
--list_file=training_set.list \
--template_strand=$template_strand \
--move_threshold=$move_threshold \
--pseudocount=$pseudocount \
--suffix_filename="${strand}_move${move_threshold}_pc${pseudocount}" \
--logtostderr

echo 'Training was successfully finished. Sampling...';

${scripts_folder}/../src/sample_move_hmm_main \
--list_file=testing_set.list \
--template_strand=$template_strand \
--trained_move_hmm="trained_move_hmm_${strand}_move${move_threshold}_pc${pseudocount}.json" \
--samples=$samples \
--logtostderr

echo 'Sampling finished!'

mv *.samples ../
