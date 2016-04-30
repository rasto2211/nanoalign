#!/bin/bash

# Run needle aligner with default options and return identity percentage.

needle --asequence $1 --bsequence $2 stdout -gapopen 10.0 -gapextend 0.5 | 
grep Identity | sed -r 's/.*\((.*)\).*/\1/'
