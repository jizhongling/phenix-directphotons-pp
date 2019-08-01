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
cp ../$0:t:r.C ../$config proc$2

cd proc$2
set pt = `echo "3 + $2/4 * 0.1" | bc`
sed -i "s/^ckin 3 .*/ckin 3 $pt/" $config
set pt = `echo "$pt + 1" | bc`
sed -i "s/^ckin 4 .*/ckin 4 $pt/" $config

root -l -b -q $0:t:r.C\($2\)

cd ..
rm -rf proc$2
