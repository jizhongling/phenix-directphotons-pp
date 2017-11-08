#!/bin/bash
# Function: Combine root files ten by ten.

#cd "$SPIN/taxi/Run13pp510MinBias/10701/data"
cd "$SPIN/taxi/Run13pp510ERT/10853/data"
rm -f total.root tmp.root

files=""
count=0

#for FILE in DirectPhotonPP-*.root ; do
#  files="${files} ${FILE}"
while read -d " " runnumber ; do
  files="${files} DirectPhotonPP-${runnumber}.root"
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
done < "$PLHF/taxi/Run13pp510ERT/runlist.txt"

if [[ -n "${files}" ]] ; then
  if [[ -f "total.root" ]] ; then
    hadd tmp.root total.root ${files}
  else
    hadd tmp.root ${files}
  fi
  mv tmp.root total.root
fi
