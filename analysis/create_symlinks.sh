#!/bin/bash

# Creates symlinks to all files in the input in the given folder.

folder=$1

mkdir $folder;
while read line;
do
    ln -s $line $folder
done
