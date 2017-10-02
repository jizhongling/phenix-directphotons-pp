#!/bin/bash
# Function: Generate "phparticlegen.txt".

if [[ ! -d "filelist" ]] ; then
  mkdir filelist
else
  echo -e "Error: Directory \"filelist\" already exists!"
  exit 1
fi

if [[ ! -d "histos" ]] ; then
  mkdir histos
fi

for i in $(seq 0 591) ; do
  fn=$(( i/20 ))
  echo "$PLHF/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/phparticlegen/phparticlegen${i}.root" >> "filelist/phparticlegen${fn}.txt"
done
