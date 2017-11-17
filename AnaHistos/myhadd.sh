#!/bin/bash
# Function: Combine root files ten by ten.

echo -en "This is for the TAXI files, do you really want to continure? [N/y] "
read -n 1 yes
echo
if [[ "${yes}" != "y" ]] ; then
  exit 0
fi

cd "$PLHF/taxi/Run13pp510ERT/12232/data"
rm -f total.root tmp.root

prename="Pi0PP-"
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

mv total.root ${prename}histo.root
