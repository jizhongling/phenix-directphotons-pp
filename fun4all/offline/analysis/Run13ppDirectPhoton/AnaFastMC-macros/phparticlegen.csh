#!/bin/tcsh
# Usage: $0 $(Initialdir) Macro.C $(Process)

setenv HOME /phenix/u/$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
  source $i
end
source $HOME/.login

cd $1
if(! -e phparticlegen) then mkdir phparticlegen
root -l -b -q $2\(20000,\"/phenix/plhf/zji/sources/offline/analysis/Run13ppDirectPhoton/AnaFastMC-macros/phparticlegen/phparticlegen$3.root\"\)
