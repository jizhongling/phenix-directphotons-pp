#!/bin/bash
# Function: Generate "simDST.txt".

if [[ ! -d "filelist" ]] ; then
  mkdir filelist
else
  echo -e "Error: Directory \"filelist\" already exists!"
  exit 1
fi

if [[ ! -d "histos" ]] ; then
  mkdir histos
fi

nfiles=10
fn=0
count=0

for i in $(seq 0 497) ; do
  if (( "${count}" >= "${nfiles}" )) ; then
    count=0
    (( fn++ ))
  fi
  echo "$SPIN/data/pisaRun13/simDST/simDST${i}.root" >> "filelist/simDST${fn}.txt"
  (( count++ ))
done
