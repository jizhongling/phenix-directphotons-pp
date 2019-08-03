#!/bin/bash
# Function: run PHPythia, PISA and pisaToDST in sequence and check output.
# Usage: $0 $(Process)

IsZombie () {
    root -l <<EOF
  {
    TFile *f = new TFile("${1}");
    if(!f || f->IsZombie()) return 1;
    TTree *t = (TTree*)f->Get("T");
    if(!t || t->GetEntries() != 10000) return 1;
    return 0;
  }
EOF
}

proc="${1}"
working="${_CONDOR_SCRATCH_DIR}/working${proc}"
output_dir="${SPIN}/data/pisaRun13/simDST-dirphoton"
logfile="${output_dir}/aa_pythiaToDST.log"
goodlist="${output_dir}/aa_goodlist.txt"
badlist="${output_dir}/aa_badlist.txt"
nProcess="4500"
scale=`echo "${nProcess}/300" | bc`
minsize="3"
maxcheck="5"

mkdir -p "${working}"
mkdir -p "${output_dir}"
cd "${working}"
echo -e "Process ${proc}: working directory ${PWD}" >> "${logfile}"

# Run PHPythia
pythia_ldir="${PLHF}/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros"
pythia_macro="anaFastMC_GenPH.C"
pythia_config="pythia.cfg"
pythia_tree="phpythia.root"
pythia_histo="${output_dir}/AnaFastMC-GenPH-histo${proc}.root"

cp "${pythia_ldir}/${pythia_macro}" .
cp "${pythia_ldir}/${pythia_config}" .

pt=`echo "3 + ${proc}/${scale} * 0.1" | bc`
sed -i "s/^ckin 3 .*/ckin 3 ${pt}/" "${pythia_config}"
pt=`echo "${pt} + 1" | bc`
sed -i "s/^ckin 4 .*/ckin 4 ${pt}/" "${pythia_config}"

echo -e "Process ${proc}: start running PHPythia" >> "${logfile}"
root -l -b -q "${pythia_macro}(${proc},${scale},\"${pythia_tree}\",\"${pythia_histo}\")"

# Test PHPythia output
for icheck in {1..10} ; do
    if [[ $(IsZombie "${pythia_tree}" | awk '{print substr($0,length($0),1)}') -eq "1" || $(ls -l --block-size=M "${pythia_tree}" | awk '{printf "%d", $5}') -lt "${minsize}" ]] ; then
	if [[ "${icheck}" -gt "${maxcheck}" ]] ; then
	    echo -e "Process ${proc}: PHPythia failed" >> "${logfile}"
	    echo -e "${proc}" >> "${badlist}"
	    exit 1
	fi
        sleep 10
	echo -e "Process ${proc}: rerun PHPythia ${icheck} times" >> "${logfile}"
        root -l -b -q "${pythia_macro}(${proc},${scale},\"${pythia_tree}\",\"${pythia_histo}\")"
    fi
done

echo -e "Process ${proc}: finish running PHPythia" >> "${logfile}"

# Run PISA
pisa_ldir="${PLHF}/data/pisaRun13"
pisa_tree="PISAEvent.root"

echo -e "phpythia 100 ${pythia_tree}" > "glogon.kumac"
echo -e "pisafile ${pisa_tree}" >> "glogon.kumac"
echo -e "ptrig 10000" >> "glogon.kumac"
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

echo -e "Process ${proc}: start running PISA" >> "${logfile}"
pisa < "pisa.input"

# Test PISA output
for icheck in {1..10} ; do
    if [[ $(IsZombie "${pisa_tree}" | awk '{print substr($0,length($0),1)}') -eq "1" || $(ls -l --block-size=M "${pisa_tree}" | awk '{printf "%d", $5}') -lt "${minsize}" ]] ; then
	if [[ "${icheck}" -gt "${maxcheck}" ]] ; then
	    echo -e "Process ${proc}: PISA failed" >> "${logfile}"
	    echo -e "${proc}" >> "${badlist}"
	    exit 1
	fi
        sleep 10
	echo -e "Process ${proc}: rerun PISA ${icheck} times" >> "${logfile}"
        pisa < "pisa.input"
    fi
done

echo -e "Process ${proc}: finish running PISA" >> "${logfile}"

# Run pisaToDST
cd "${pisa_ldir}"
pisa_macro="pisaToDST.C"
pisa_tree="${working}/${pisa_tree}"
dst_tree="${output_dir}/simDST${proc}.root"

echo -e "Process ${proc}: start running pisaToDST" >> "${logfile}"
root -l -b -q "${pisa_macro}(0,\"${pisa_tree}\",\"${dst_tree}\")"

# Test pisaToDST output
for icheck in {1..10} ; do
    if [[ $(IsZombie "${dst_tree}" | awk '{print substr($0,length($0),1)}') -eq "1" || $(ls -l --block-size=M "${dst_tree}" | awk '{printf "%d", $5}') -lt "${minsize}" ]] ; then
	if [[ "${icheck}" -gt "${maxcheck}" ]] ; then
	    echo -e "Process ${proc}: pisaToDST failed" >> "${logfile}"
	    echo -e "${proc}" >> "${badlist}"
	    exit 1
	fi
        sleep 10
	echo -e "Process ${proc}: rerun pisaToDST ${icheck} times" >> "${logfile}"
        root -l -b -q "${pisa_macro}(0,\"${pisa_tree}\",\"${dst_tree}\")"
    fi
done

echo -e "Process ${proc}: finish running pisaToDST" >> "${logfile}"
echo -e "${proc}" >> "${goodlist}"
