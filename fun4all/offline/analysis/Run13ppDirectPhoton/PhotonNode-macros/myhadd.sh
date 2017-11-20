#!/bin/bash
# Function: Combine root files ten by ten.

cd "histos-ERT"
rm -f total.root tmp.root

prename="PhotonNode-"
files=""
count=0

#for FILE in ${prename}*.root ; do
while read -d " " runnumber ; do
  files="${files} ${prename}${runnumber}.root"
  (( count++ ))
  if (( "${count}" > "9" )) ; then
    if [[ -f "total.root" ]] ; then
      hadd tmp.root total.root ${files}
    else
      hadd tmp.root ${files}
    fi
    mv tmp.root total.root
    files=""
    count=0
  fi
done < "/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt"

if [[ -n "${files}" ]] ; then
  if [[ -f "total.root" ]] ; then
    hadd tmp.root total.root ${files}
  else
    hadd tmp.root ${files}
  fi
  mv tmp.root total.root
fi
