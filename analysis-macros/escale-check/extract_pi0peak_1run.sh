#!/bin/zsh

# single run test example
RUNINDEX=1
RUNNUMBER=387565

HISTFILE=/gpfs/mnt/gpfs02/phenix/spin3/nfeege/taxi/Run13pp510ERT/8401/data/DirectPhotonPP-387565.root
HISTNAME=inv_mass_2photon_pi0calib_sector

WRITEPLOTS=1

root -l extract_pi0peak.C\($RUNINDEX,$RUNNUMBER,\"$HISTFILE\",\"$HISTNAME\",$WRITEPLOTS\)

exit
