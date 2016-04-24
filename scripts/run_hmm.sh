#!/bin/bash

scripts_folder=~/testing_scripts/nanopore-read-align/scripts
reads_folder=~/ecoli
working_folder=~/2d_reads_move2_training90
strand=template
template_strand=true
training_set_percent=0.9
move_threshold=2
pseudocount=1
samples=256

mkdir -p $working_folder
cd $working_folder;
mkdir -p hmm;
cd hmm;

echo 'Filtering reads...';

python3 ${scripts_folder}/list_2d_reads_with_move_threshold.py \
--reads_folder=$reads_folder \
--move_threshold=$move_threshold \
--strand=$strand > 2d_reads_${strand}_move${move_threshold}.list

echo 'Finished filtering. Creating symlinks...'

while read line;
do
    ln -s $line ../
done < 2d_reads_${strand}_move${move_threshold}.list;

echo 'Finished creating symlinks. Splitting data...';

python3 ${scripts_folder}/split_data.py \
--file_with_list="2d_reads_${strand}_move${move_threshold}.list" \
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
