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

set NFiles = 2
@ START = $2 * $NFiles
@ END = ( $2 + 1 ) * $NFiles - 1

setenv INPUT
foreach proc ( `seq $START $END` )
  foreach i ( `seq 1 2` )
    setenv GZ $SPIN/data/powheg/pwgevents$proc-`printf "%04d" $i`.lhe.gz
    setenv LHE $SPIN/data/powheg/proc$proc/pwgevents-`printf "%04d" $i`.lhe
    if ( -f $GZ ) then
      if ( `ls -l --block-size=M $GZ | awk '{printf "%d", $5}'` > 3 ) then
        setenv INPUT "$INPUT $GZ"
      endif
    else if ( -f $LHE ) then
      if ( `ls -l --block-size=M $LHE | awk '{printf "%d", $5}'` > 3 ) then
        setenv INPUT "$INPUT $LHE"
      endif
    endif
  end
end

cd $1
./anaLHEF anaLHEF.cmnd histos/AnaPowheg-histo$2.root $INPUT
./anaLHEF anaLHEFNoMPI.cmnd histos/AnaPowhegNoMPI-histo$2.root $INPUT
