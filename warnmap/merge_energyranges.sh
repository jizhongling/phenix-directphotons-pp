#!/bin/zsh

FILE1="warnmap-data/WarnmapData_Run13pp510MinBias.root"
FILE2="warnmap-data/WarnmapData_Run13pp510MinBias_ERange.root"
root -b -q merge_energyranges.C\(\"${FILE1}\",\"${FILE2}\"\)

FILE1="warnmap-data/WarnmapData_Run13pp510MinBias_RunRange1.root"
FILE2="warnmap-data/WarnmapData_Run13pp510MinBias_RunRange1_ERange.root"
root -b -q merge_energyranges.C\(\"${FILE1}\",\"${FILE2}\"\)

FILE1="warnmap-data/WarnmapData_Run13pp510MinBias_RunRange2.root"
FILE2="warnmap-data/WarnmapData_Run13pp510MinBias_RunRange2_ERange.root"
root -b -q merge_energyranges.C\(\"${FILE1}\",\"${FILE2}\"\)

FILE1="warnmap-data/WarnmapData_Run13pp510MinBias_RunRange3.root"
FILE2="warnmap-data/WarnmapData_Run13pp510MinBias_RunRange3_ERange.root"
root -b -q merge_energyranges.C\(\"${FILE1}\",\"${FILE2}\"\)

FILE1="warnmap-data/WarnmapData_Run13pp510MinBias_RunRange4.root"
FILE2="warnmap-data/WarnmapData_Run13pp510MinBias_RunRange4_ERange.root"
root -b -q merge_energyranges.C\(\"${FILE1}\",\"${FILE2}\"\)

FILE1="warnmap-data/WarnmapData_Run13pp510MinBias_RunRange5.root"
FILE2="warnmap-data/WarnmapData_Run13pp510MinBias_RunRange5_ERange.root"
root -b -q merge_energyranges.C\(\"${FILE1}\",\"${FILE2}\"\)

