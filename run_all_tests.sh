#!/bin/bash

cd tests;
for test_file in *_test;
do
  ./$test_file | grep 'FAILED'
done
