#!/bin/bash
# Function: Concatenate txt files.

for ipt in $(seq 0 29) ; do
  cat ALL/Process*-pT${ipt}.txt > ALL/pT${ipt}.txt
done

rm ALL/Process*.txt
