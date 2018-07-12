#!/bin/tcsh
# Usage: $0 $(Initialdir) $(Process)

setenv HOME /phenix/u/$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
  source $i
end
source $HOME/.login

set config="pythia.cfg"

cd $1/histos
mkdir proc$2
cp ../phpythia.C ../$config proc$2

cd proc$2
set pt = `echo "5 + $2/5 * 0.1" | bc`
sed -i "s/ckin 3 .*/ckin 3 $pt/" $config

root -l -b -q phpythia.C\(500000,\"../AnaFastMC-PH-histo$2.root\"\)

cd ..
rm -rf proc$2
