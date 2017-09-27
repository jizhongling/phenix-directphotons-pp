#!/bin/tcsh
# Usage: $0 $(Initialdir) $(Process)

setenv HOME /phenix/u/$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
  source $i
end
source $HOME/.login

cd $1
root -l -b -q $0:t:r.C\(\"filelist/phparticlegen$2.txt\",\"histos$2.root\"\)