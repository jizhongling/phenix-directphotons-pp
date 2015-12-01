#!/bin/zsh

root -b -q plot_warnmap2D.C\(\"warnmap-final/Warnmap_Run13pp510MinBias_Final.txt\",true\)
root -b -q plot_warnmap2D.C\(\"warnmap-final/Warnmap_Run13pp510ERT_Final.txt\",true\)
root -b -q plot_warnmap2D.C\(\"warnmap-final/Warnmap_Run9pp500MinBias_Final.txt\",true\)
root -b -q plot_warnmap2D.C\(\"warnmap-final/Warnmap_Run9pp500ERT_Final.txt\",true\)
root -b -q plot_warnmap2D.C\(\"warnmap-final/Warnmap_Run9pp200MinBias_Final.txt\",true\)
root -b -q plot_warnmap2D.C\(\"warnmap-final/Warnmap_Run9pp200ERT_Final.txt\",true\)

root -b -q plot_warnmap2D.C\(\"warnmap-output/Warnmap_Run13pp510MinBias_erange_0_nsigma10_niter10.txt\",true\)
root -b -q plot_warnmap2D.C\(\"warnmap-output/Warnmap_Run13pp510MinBias_erange_4_nsigma10_niter10.txt\",true\)

root -b -q plot_warnmap2D.C\(\"warnmap-output/Warnmap_Run13pp510MinBias_RunRange5_ERange_erange_0_to_4_nsigma10_niter10.txt\",true\)
root -b -q plot_warnmap2D.C\(\"warnmap-output/Warnmap_Run13pp510MinBias_RunRange4_ERange_erange_0_to_4_nsigma10_niter10.txt\",true\)

