#!/bin/bash
# Function: combine root files N by N.

set -o errexit
set -o nounset
set -o pipefail

cd "histos-TAXI"
rm -f total.root tmp.root

prename="PhotonHistos-"
files=""
count=0

for FILE in ${prename}*.root ; do
    files="${files} ${FILE}"
    (( ++count ))
    if (( "${count}" >= "10" )) ; then
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

#mv total.root "../${prename}Sasha.root"
#mv total.root "../${prename}DC3sigma.root"
mv total.root "../${prename}Inseok-tightcut.root"
