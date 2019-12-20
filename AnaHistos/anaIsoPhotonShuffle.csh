#!/bin/tcsh
# Usage: $0 $(Initialdir) $(Process)

setenv HOME /phenix/u/$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
  source $i
end
source $HOME/.login

set TAXI = 15763
set nProcess = 20
set proc = `echo "$2 * $nProcess" | bc`
set proc_end = `echo "($2+1) * $nProcess" | bc`

cd $1
while ( $proc < $proc_end )
  ./anaIsoPhotonShuffle $TAXI $proc
  @ proc++
end
