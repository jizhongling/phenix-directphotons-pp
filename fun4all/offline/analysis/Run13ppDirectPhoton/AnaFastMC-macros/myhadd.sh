#!/bin/bash
# Function: Combine root files N by N.

set -o errexit
set -o nounset
set -o pipefail

cd "histos"
rm -f total.root tmp.root

prename="AnaFastMC-PH-"
files=""
count=0

for FILE in ${prename}*.root ; do
    files="${files} ${FILE}"
    (( ++count ))
    if (( "${count}" >= "100" )) ; then
        if [[ -f "total.root" ]] ; then
    	    hadd tmp.root total.root ${files}
        else
    	    hadd tmp.root ${files}
        fi
        mv tmp.root total.root
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
    mv tmp.root total.root
fi

mv total.root ../${prename}histo.root
