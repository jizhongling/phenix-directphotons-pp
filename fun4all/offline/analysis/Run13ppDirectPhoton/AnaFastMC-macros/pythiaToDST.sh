#!/bin/bash
# Function: run PHPythia, PISA and pisaToDST in sequence and check output.
# Usage: $0 $(Process)

IsZombie () {
    root -l <<EOF
  {
    TFile f("${1}");
    if(f.IsZombie()) return 1;
    else return 0;
  }
EOF
}

proc="${1}"
working="${_CONDOR_SCRATCH_DIR}/working${proc}"
nProcess=12000
scale=`echo "${nProcess}/300" | bc`
output_dir="${SPIN}/data/pisaRun13/simDST-dirphoton-0"

mkdir -p "${working}"
mkdir -p "${output_dir}"
cd "${working}"

# Run PHPythia
pythia_ldir="${PLHF}/github/phenix-directphotons-pp/fun4all/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros"
pythia_macro="anaFastMC_GenPH.C"
pythia_config="pythia.cfg"
pythia_tree="phpythia.root"
pythia_histo="${output_dir}/AnaFastMC-GenPH-histo${proc}.root"
pythia_minMB=8

cp "${pythia_ldir}/${pythia_macro}" .
cp "${pythia_ldir}/${pythia_config}" .

pt=`echo "3 + ${proc}/${scale} * 0.1" | bc`
sed -i "s/^ckin 3 .*/ckin 3 ${pt}/" "${pythia_config}"
pt=`echo "${pt} + 1" | bc`
sed -i "s/^ckin 4 .*/ckin 4 ${pt}/" "${pythia_config}"

root -l -b -q "${pythia_macro}(${proc},\"${pythia_tree}\",\"${pythia_histo}\")"

# Test PHPythia output
for icheck in {1..2} ; do
    if [[ $(IsZombie "${pythia_tree}" | awk '{print substr($0,length($0),1)}') = 1 || $(ls -l --block-size=M "${pythia_tree}" | awk '{printf "%d", $5}') < ${pythia_minMB} ]] ; then
        sleep 10
        root -l -b -q "${pythia_macro}(${proc},\"${pythia_tree}\",\"${pythia_histo}\")"
    fi
done

# Run PISA
pisa_ldir="${PLHF}/data/pisaRun13"
pisa_tree="PISAEvent.root"
pisa_minKB=8

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

pisa < "pisa.input"

# Test PISA output
for icheck in {1..2} ; do
    if [[ $(IsZombie "${pisa_tree}" | awk '{print substr($0,length($0),1)}') = 1 || $(ls -l --block-size=K "${pisa_tree}" | awk '{printf "%d", $5}') < ${pisa_minKB} ]] ; then
        sleep 10
        pisa < "pisa.input"
    fi
done

# Run pisaToDST
cd "${pisa_ldir}"
pisa_macro="pisaToDST.C"
pisa_tree="${working}/PISAEvent.root"
dst_tree="${output_dir}/simDST${proc}.root"
dst_minKB=8

root -l -b -q "${pisa_macro}(0,\"${pisa_tree}\",\"${dst_tree}\")"

# Test pisaToDST output
for icheck in {1..2} ; do
    if [[ $(IsZombie "${dst_tree}" | awk '{print substr($0,length($0),1)}') = 1 || $(ls -l --block-size=K "${dst_tree}" | awk '{printf "%d", $5}') < ${dst_minKB} ]] ; then
        sleep 10
        root -l -b -q "${pisa_macro}(0,\"${pisa_tree}\",\"${dst_tree}\")"
    fi
done

rm -rf "${working}"
