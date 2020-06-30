#!/bin/tcsh
# Usage: $0 $(Initialdir) $(Process)

setenv HOME /phenix/u/$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
  source $i
end
source $HOME/bin/phenix_setup_64.csh -n new

setenv MYINSTALL /phenix/plhf/zji/install
setenv PATH $HOME/bin:$MYINSTALL/bin:$PATH
setenv LD_LIBRARY_PATH $MYINSTALL/lib64:$MYINSTALL/lib:$LD_LIBRARY_PATH
setenv TSEARCHPATH "$MYINSTALL":$TSEARCHPATH

setenv PLHF /phenix/plhf/zji
setenv SPIN /phenix/spin/phnxsp01/zji
setenv SCRATCH /phenix/scratch/zji

cd $1
foreach i (`seq 1 2`)
  ./anaLHEF histos/powheg-histo$2-`printf "%04d" $i`.root $SPIN/data/powheg/pwgevents$2-`printf "%04d" $i`.lhe.gz &
end
wait

echo Finished
