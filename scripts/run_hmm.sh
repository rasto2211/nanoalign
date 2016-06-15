#!/bin/bash

# Takes file with list of fast5 files and path of working folder.
# All files produced by HMM will be placed to working folder.
# Files with list of files used for training/testing and trained HMM
# are placed into separate subfolder $working_folder/hmm.

# Constants declaration.

training_set_percent=0.7 # Split data to training and testing.
move_threshold=2
pseudocount=1
samples=250
# Folder with all the scripts.
scripts_folder=~/nanopore-read-align/scripts

input_reads_list=$1
working_folder=$2
strand=template
template_strand=true # Template vs complement

# End of constants declarations.

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
