#!/bin/bash
# Function: Concatenate txt files.

for itype in $(seq 0 2) ; do
  for ipart in $(seq 0 2) ; do
    for icr in $(seq 0 1) ; do
      for ipt in $(seq 0 29) ; do
        cat ALL/Process*-type${itype}-part${ipart}-crossing${icr}-pT${ipt}.txt > ALL/type${itype}-part${ipart}-crossing${icr}-pT${ipt}.txt
      done
    done
  done
done

rm ALL/Process*.txt
