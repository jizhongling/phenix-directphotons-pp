#!/bin/zsh

checkfiles=(
'warnmap-output/Checkplots_Run13pp510MinBias_ybins3to4_nsigma10_niter10.root'
'warnmap-output/Checkplots_Run13pp510MinBias_ybins5to6_nsigma10_niter10.root'
'warnmap-output/Checkplots_Run13pp510MinBias_ybins7to9_nsigma10_niter10.root'
'warnmap-output/Checkplots_Run13pp510MinBias_ybins10to12_nsigma10_niter10.root'
'warnmap-output/Checkplots_Run13pp510MinBias_ybins13to25_nsigma10_niter10.root'
'warnmap-output/Checkplots_Run13pp510MinBias_ybins3to25_nsigma10_niter10.root' )

for file in $checkfiles; do

    echo "Next up: $file"

    root -b -q plot_hitrate_threshold.C\(\"$file\",1\)

done



#root -b -q plot_hitrate_threshold.C\(\"warnmap-output/Checkplots_Run13pp510MinBias_ybins3to4_nsigma10_niter10.root\",1\)