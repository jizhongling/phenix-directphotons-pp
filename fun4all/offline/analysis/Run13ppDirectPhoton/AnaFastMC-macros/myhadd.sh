#!/bin/bash
# Function: combine good root files N by N.

set -o errexit
set -o nounset
set -o pipefail

output_dir="${PWD}"
cd "$SPIN/data/pythiaToHisto"
rm -f pythia_total.root pythia_tmp.root
rm -f dst_total.root dst_tmp.root

pythia_prename="AnaFastMC-GenPH-histo"
dst_prename="HadronResponse-histo"
pythia_files=""
dst_files=""
count=0

while read -d " " runnumber ; do
    pythia_files="${pythia_files} ${pythia_prename}${runnumber}.root"
    dst_files="${dst_files} ${dst_prename}${runnumber}.root"
    (( ++count ))
    if (( "${count}" >= "900" )) ; then
	if [[ -f "pythia_total.root" ]] ; then
	    hadd pythia_tmp.root pythia_total.root ${pythia_files}
	else
	    hadd pythia_tmp.root ${pythia_files}
	fi
	if [[ -f "dst_total.root" ]] ; then
	    hadd dst_tmp.root dst_total.root ${dst_files}
	else
	    hadd dst_tmp.root ${dst_files}
	fi
	mv pythia_tmp.root pythia_total.root
	mv dst_tmp.root dst_total.root
	pythia_files=""
	dst_files=""
	count=0
    fi
done < "aa_goodlist.txt"

if [[ -n "${pythia_files}" || -n "${dst_files}" ]] ; then
    if [[ -f "pythia_total.root" ]] ; then
	hadd pythia_tmp.root pythia_total.root ${pythia_files}
    else
	hadd pythia_tmp.root ${pythia_files}
    fi
    if [[ -f "dst_total.root" ]] ; then
	hadd dst_tmp.root dst_total.root ${dst_files}
    else
	hadd dst_tmp.root ${dst_files}
    fi
    mv pythia_tmp.root pythia_total.root
    mv dst_tmp.root dst_total.root
fi

mv pythia_total.root ${output_dir}/${pythia_prename}.root
mv dst_total.root ${output_dir}/${dst_prename}.root
