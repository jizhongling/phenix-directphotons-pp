#!/bin/bash
# Function: run PHPythia, PISA, pisaToDST and anaDST in sequence and check output.
# Usage: $0 $(Process)

# General variables
proc="${1}"
working="${_CONDOR_SCRATCH_DIR}/working${proc}"
output_dir="${SPIN}/data/pythiaToHisto"
logfile="${output_dir}/aa_pythiaToHisto.log"
goodlist="${output_dir}/aa_goodlist.txt"
badlist="${output_dir}/aa_badlist.txt"
nProcess="12000"
#nProcess="1200"
scale=`echo "${nProcess}/300" | bc`
nevents="10000"
#nevents="5000"
maxcheck="4"
minsize="1"

# PHPythia variables
pythia_ldir="${PLHF}/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros"
file_checker="${pythia_ldir}/IsGoodFile"
pythia_macro="anaFastMC_GenPH.C"
#pythia_macro="phparticlegen.C"
pythia_config="pythia.cfg"
pythia_xsec="${output_dir}/phpy_xsec${proc}.root"
pythia_tree="phpythia${proc}.root"
pythia_histo="AnaFastMC-GenPH-histo${proc}.root"

# PISA variables
pisa_ldir="${PLHF}/data/pisaRun13"
pisa_tree="PISAEvent${proc}.root"

# pisaToDST variables
pisa_macro="pisaToDST.C"
dst_tree="simDST${proc}.root"

# anaDST variables
dst_macro="anaDST.C"
dst_histo="HadronResponse-histo${proc}.root"

# Exit if already done
if [[ -f "${output_dir}/${dst_histo}" ]] ; then
    exit 0
fi

# Prepare to run
mkdir -p "${working}"
mkdir -p "${output_dir}"
cd "${working}"
cp "${pythia_xsec}" "phpy_xsec.root"
cp "${output_dir}/${pythia_histo}" .
cp "${output_dir}/${pythia_tree}" .
cp "${output_dir}/${pisa_tree}" .
cp "${output_dir}/${dst_tree}" .

# Run PHPythia
cd "${working}"
cp "${pythia_ldir}/${pythia_macro}" .
cp "${pythia_ldir}/${pythia_config}" .

pt=`echo "3 + ${proc}/${scale} * 0.1" | bc`
sed -i "s/^ckin 3 .*/ckin 3 ${pt}/" "${pythia_config}"
pt=`echo "${pt} + 1" | bc`
sed -i "s/^ckin 4 .*/ckin 4 ${pt}/" "${pythia_config}"

for (( icheck = 1; icheck <= ${maxcheck}; icheck++ )) ; do
    if [[ -f "${pisa_tree}" || -f "${dst_tree}" ]] ; then
        break
    fi
    "${file_checker}" "histo" "${pythia_histo}" "${nevents}"
    status_histo="${?}"
    #status_histo="0"
    "${file_checker}" "tree" "${pythia_tree}" "${nevents}"
    status_tree="${?}"
    if [[ "${status_histo}" -ne "0" || "${status_tree}" -ne "0" || $(ls -l --block-size=M "${pythia_tree}" | awk '{printf "%d", $5}') -lt "${minsize}" ]] ; then
	if [[ "${icheck}" -eq "${maxcheck}" ]] ; then
	    echo -e "Process ${proc}: PHPythia failed" >> "${logfile}"
	    echo -en "${proc} " >> "${badlist}"
	    exit 1
	fi
        sleep 10
	echo -e "Process ${proc}: PHPythia runs ${icheck} times" >> "${logfile}"
        root -l -b -q "${pythia_macro}(${proc},${scale},\"${pythia_tree}\",\"${pythia_histo}\")"
        #root -l -b -q "${pythia_macro}(${proc},${nevents},\"${pythia_tree}\")"
    else
        rm -f "${output_dir}/${pythia_tree}"
        echo -e "Process ${proc}: PHPythia finished running" >> "${logfile}"
        break
    fi
done

# Run PISA
cd "${working}"

echo -e "phpythia 100 ${pythia_tree}" > "glogon.kumac"
echo -e "pisafile ${pisa_tree}" >> "glogon.kumac"
echo -e "ptrig ${nevents}" >> "glogon.kumac"
echo -e "exit" >> "glogon.kumac"

cp "${pisa_ldir}/pisa.kumac" "pisa.kumac"
let "rnd = ${proc} % 215 + 1"
sed -i "s/RNDM 001 0/RNDM ${rnd} 0/" "pisa.kumac"

ln -s "${pisa_ldir}/pisa.input"
ln -s "${pisa_ldir}/gffgo.dat"
ln -s "${pisa_ldir}/phnx.par"
ln -s "${pisa_ldir}/event.par"
ln -s "${pisa_ldir}/flukaaf.dat"
ln -s "${pisa_ldir}/xsneut95.dat"
ln -s "${pisa_ldir}/Sim3D++.root"

for (( icheck = 1; icheck <= ${maxcheck}; icheck++ )) ; do
    if [[ -f "${dst_tree}" ]] ; then
        break
    fi
    "${file_checker}" "tree" "${pisa_tree}" "${nevents}"
    status_tree="${?}"
    if [[ "${status_tree}" -ne "0" || $(ls -l --block-size=M "${pisa_tree}" | awk '{printf "%d", $5}') -lt "${minsize}" ]] ; then
	if [[ "${icheck}" -eq "${maxcheck}" ]] ; then
	    cp "phpy_xsec.root" "${pythia_xsec}"
	    cp "${pythia_histo}" "${output_dir}"
            cp "${pythia_tree}" "${output_dir}"
	    echo -e "Process ${proc}: PISA failed" >> "${logfile}"
	    echo -en "${proc} " >> "${badlist}"
	    exit 1
	fi
        sleep 10
	echo -e "Process ${proc}: PISA runs ${icheck} times" >> "${logfile}"
        pisa < "pisa.input"
    else
        rm -f "${output_dir}/${pythia_tree}"
        rm -f "${output_dir}/${pisa_tree}"
        echo -e "Process ${proc}: PISA finished running" >> "${logfile}"
        break
    fi
done

# Run pisaToDST
cd "${pisa_ldir}"

for (( icheck = 1; icheck <= ${maxcheck}; icheck++ )) ; do
    "${file_checker}" "tree" "${working}/${dst_tree}" "${nevents}"
    status_tree="${?}"
    if [[ "${status_tree}" -ne "0" || $(ls -l --block-size=M "${working}/${dst_tree}" | awk '{printf "%d", $5}') -lt "${minsize}" ]] ; then
	if [[ "${icheck}" -eq "${maxcheck}" ]] ; then
	    cp "${working}/phpy_xsec.root" "${pythia_xsec}"
	    cp "${working}/${pythia_histo}" "${output_dir}"
            cp "${working}/${pisa_tree}" "${output_dir}"
            rm -f "${output_dir}/${pythia_tree}"
	    echo -e "Process ${proc}: pisaToDST failed" >> "${logfile}"
	    echo -en "${proc} " >> "${badlist}"
	    exit 1
	fi
	sleep 10
	echo -e "Process ${proc}: pisaToDST runs ${icheck} times" >> "${logfile}"
	root -l -b -q "${pisa_macro}(0,\"${working}/${pisa_tree}\",\"${working}/${dst_tree}\")"
    else
        rm -f "${output_dir}/${pythia_tree}"
        rm -f "${output_dir}/${pisa_tree}"
        rm -f "${output_dir}/${dst_tree}"
        echo -e "Process ${proc}: pisaToDST finished running" >> "${logfile}"
	break
    fi
done

# Run anaDST
cd "${working}"
cp "${pythia_ldir}/${dst_macro}" .

for (( icheck = 1; icheck <= ${maxcheck}; icheck++ )) ; do
    "${file_checker}" "histo" "${dst_histo}" "${nevents}"
    status_histo="${?}"
    if [[ "${status_histo}" -ne "0" ]] ; then
	if [[ "${icheck}" -eq "${maxcheck}" ]] ; then
	    cp "phpy_xsec.root" "${pythia_xsec}"
	    cp "${pythia_histo}" "${output_dir}"
            cp "${dst_tree}" "${output_dir}"
            rm -f "${output_dir}/${pythia_tree}"
            rm -f "${output_dir}/${pisa_tree}"
	    echo -e "Process ${proc}: anaDST failed" >> "${logfile}"
	    echo -en "${proc} " >> "${badlist}"
	    exit 1
	fi
        sleep 10
	echo -e "Process ${proc}: anaDST runs ${icheck} times" >> "${logfile}"
        root -l -b -q "${dst_macro}(${proc},${scale},\"${dst_tree}\",\"${dst_histo}\")"
    else
        rm -f "${pythia_xsec}"
        rm -f "${output_dir}/${pythia_tree}"
        rm -f "${output_dir}/${pisa_tree}"
        rm -f "${output_dir}/${dst_tree}"
        echo -e "Process ${proc}: anaDST finished running" >> "${logfile}"
        break
    fi
done

# Store histograms
cd "${working}"
cp "${pythia_histo}" "${output_dir}"
cp "${dst_histo}" "${output_dir}"
echo -en "${proc} " >> "${goodlist}"
