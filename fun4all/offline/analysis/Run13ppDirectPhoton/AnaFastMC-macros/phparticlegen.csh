#!/bin/tcsh
# Usage: $0 $(Initialdir) $(Process)

setenv HOME /phenix/u/$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
  source $i
end
source $HOME/.login

cd $1
root -l -b -q $0:t:r.C\(20000,\"/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/phparticlegen/phparticlegen$2.root\"\)
