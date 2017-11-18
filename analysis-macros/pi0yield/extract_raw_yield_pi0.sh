#!/bin/zsh

DATASETS=(
    Run13pp510ERT
)

INPUT_DIR_BASE="/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/nfeege/taxi/"
OUTPUT_DIR="./"

RUNQALIST="../runqa/Run13pp510_RunQuality.txt"
#RUNQALIST="RunQA_Run13_Test.txt"

WRITEPLOTS=0


#### START DEBUG
#RUNINDEX=1
#RUNNUMBER=387431
#HISTFILE="taxi/Run13pp510MinBias/12054/data/DirectPhotonPP-387431.root"
##HISTNAME="h3_inv_mass_pi0calib_notof"
#HISTNAME="h3_inv_mass_pi0calib"
#
#root -l extract_raw_yield_pi0.C\($RUNINDEX,$RUNNUMBER,\"$HISTFILE\",\"$HISTNAME\",$WRITEPLOTS\)
#### END DEBUG


for DATASET in $DATASETS; do

    echo "Processing data set $DATASET ..."

    OUTPUT_FILE="${OUTPUT_DIR}/raw_pi0_yield_${DATASET}.txt"

    if [[ -e $OUTPUT_FILE ]]; then
	rm $OUTPUT_FILE
	touch $OUTPUT_FILE
    fi

    INPUT_DIR="$INPUT_DIR_BASE/$DATASET"

    [[ ! -d $INPUT_DIR ]] && echo "ERROR: Directory $INPUT_DIR does not exist. EXIT." && exit

    TAXI_DATA_DIR=`awk '$2=="Run_DirectPhotonPP" {print $1}' $INPUT_DIR/taxilist.txt`

    INPUT_DIR="$INPUT_DIR/$TAXI_DATA_DIR/data/"

    echo "Checking pi0 peaks from files in input directory $INPUT_DIR ... "

    FILELIST=("${(@f)$(ls $INPUT_DIR/DirectPhotonPP*.root)}")

    RUNINDEX=0
    for FILE in $FILELIST; do

	RUNNUMBER=$(echo $FILE | cut -d\- -f2 | cut -d. -f1)

	# check run QA
	RUNQA=$(cat $RUNQALIST | grep "$RUNNUMBER" | cut -d\  -f2)

	echo "Run $RUNNUMBER -> QA value $RUNQA"

	if [[ $RUNQA == 0 ]]; then
	    RUNINDEX=$(($RUNINDEX+1))

	    HISTFILE=$FILE
	    HISTNAME=hn_pi0

	    root -b -q extract_raw_yield_pi0.C\($RUNINDEX,$RUNNUMBER,\"$HISTFILE\",\"$HISTNAME\",$WRITEPLOTS\) | grep "^Raw pi0 yield count" >> $OUTPUT_FILE
	    #root -b -q extract_raw_yield_pi0.C\($RUNINDEX,$RUNNUMBER,\"$HISTFILE\",\"$HISTNAME\",$WRITEPLOTS\)
	fi

    done

done

exit
