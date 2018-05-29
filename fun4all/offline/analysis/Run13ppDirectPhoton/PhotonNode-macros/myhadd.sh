#!/bin/bash
# Function: Combine 20 root files.
# Usage: $0 $(Process)

set -o errexit
set -o nounset
set -o pipefail

cd "histos-ERT"

prename="PhotonNode-"
files=""
count=0
start=$(( $1 * 20 ))
end=$(( ($1+1) * 20 ))

#for FILE in ${prename}*.root ; do
while read -d " " runnumber ; do
  if (( "${count}" < "${start}" )) ; then
    (( ++count ))
    continue
  elif (( "${count}" >= "${end}" )) ; then
    break
  fi
  files="${files} ${prename}${runnumber}.root"
  (( ++count ))
done < "/phenix/plhf/zji/taxi/Run13pp510MinBias/runlist.txt"

hadd -f "tmp-$1.root" ${files}
