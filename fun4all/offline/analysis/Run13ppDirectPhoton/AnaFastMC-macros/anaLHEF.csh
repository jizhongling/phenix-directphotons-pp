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

@ START = $2 * 10
@ END = $2 * 10 + 9

cd $1
foreach i ( `seq 1 2` )

  setenv INPUT
  foreach proc ( `seq $START $END` )
    setenv FILE $SPIN/data/powheg/pwgevents$proc-`printf "%04d" $i`.lhe.gz
    if ( -f $FILE ) then
      if ( `ls -l --block-size=M $FILE | awk '{printf "%d", $5}'` > 50 ) then
        setenv INPUT "$INPUT $FILE"
      endif
    endif
  end

  ./anaLHEF histos/AnaPowheg-histo$2-`printf "%04d" $i`.root $INPUT &
end
wait

echo Finished
