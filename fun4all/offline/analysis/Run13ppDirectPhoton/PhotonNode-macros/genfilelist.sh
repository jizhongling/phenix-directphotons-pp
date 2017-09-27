#!/bin/bash
# Function: Generate "taxifiles.txt".

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

while read -d " " runnumber ; do
  if (( "${count}" >= "${nfiles}" )) ; then
    count=0
    (( fn++ ))
  fi
  #echo "$SPIN/taxi/Run13pp510ERT/11465/data/DirectPhotonPP_PhotonNode-${runnumber}.root" >> "filelist/taxifiles${fn}.txt"
  echo "$PLHF/taxi/Run13pp510MinBias/11343/data/DirectPhotonPP_PhotonNode-${runnumber}.root" >> "filelist/taxifiles${fn}.txt"
  (( count++ ))
done < "$PLHF/taxi/Run13pp510ERT/runlist.txt"
