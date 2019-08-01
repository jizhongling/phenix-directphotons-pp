#!/bin/bash
# Function: Combine $NFiles root files.
# Usage: $0 $(Process)

set -o errexit
set -o nounset
set -o pipefail

outdir="$(pwd)/histos-TAXI"
cd "$SPIN/taxi/Run13pp510ERT/15226/data"
#cd "$SPIN/taxi/Run13pp510MinBias/14532/data"

prename="PhotonHistos-"
NFiles=10
files=""
count=0
start=$(( $1 * ${NFiles} ))
end=$(( ($1+1) * ${NFiles} ))

while read -d " " runnumber ; do
  if (( "${count}" < "${start}" )) ; then
    (( ++count ))
    continue
  elif (( "${count}" >= "${end}" )) ; then
    break
  fi
#if [[ -e "${prename}${runnumber}.root" ]] ; then
  files="${files} ${prename}${runnumber}.root"
#fi
  (( ++count ))
done < "$PLHF/taxi/Run13pp510MinBias/runlist-DC3sigma.txt"
#done < "$PLHF/taxi/Run13pp510MinBias/runlist-Sasha.txt"

hadd -f "${outdir}/${prename}$1.root" ${files}
