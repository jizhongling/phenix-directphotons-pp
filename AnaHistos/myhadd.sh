#!/bin/bash
# Function: Combine root files ten by ten.

cd "$PLHF/taxi/Run13pp510ERT/8511/data"
rm -f total.root tmp.root

files=""
count=0

for FILE in DirectPhotonPP-*.root ; do
  if (( "${count}" < "10" )) ; then
    files="${files} ${FILE}"
    count=$(( count + 1 ))
  else
    if [[ -f "total.root" ]] ; then
      hadd tmp.root total.root ${files}
    else
      hadd tmp.root ${files}
    fi
    mv -f tmp.root total.root
    files=""
    count=0
  fi
done

if [[ -n "${files}" ]] ; then
  if [[ -f "total.root" ]] ; then
    hadd tmp.root total.root ${files}
  else
    hadd tmp.root ${files}
  fi
  mv -f tmp.root total.root
fi
