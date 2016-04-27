#!/bin/bash

cd tests;
tests=*_test;
num_tests=`echo $tests | wc -w`

tput bold
printf "Running %d tests...\n\n" $num_tests
tput sgr0

for test_file in $tests;
do
  out=`./$test_file | grep 'FAILED'`

  if [[ -z $out ]]; then
    tput setaf 2 # green color
    printf "$test_file OK\n"
  else
    tput setaf 1 # red color
    printf "$test_file FAIL\n"
    echo "$out" 
  fi

done
