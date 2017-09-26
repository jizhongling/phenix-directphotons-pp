#!/bin/zsh

DATASETS=(
    Run13pp510ERT
)

INPUT_DIR_BASE="/phenix/spin3/nfeege/taxi/"
OUTPUT_DIR="./"

RUNQALIST="../runqa/Run13pp510_RunQuality.txt"

WRITEPLOTS=0

for DATASET in $DATASETS; do

    echo "Processing data set $DATASET ..."

    OUTPUT_FILE="${OUTPUT_DIR}/pi0peak_fit_${DATASET}.txt"
    OUTPUT_FILE_RAW="${OUTPUT_DIR}/pi0peak_fit_raw_${DATASET}.txt"

    if [[ -e $OUTPUT_FILE ]]; then
	rm $OUTPUT_FILE
	touch $OUTPUT_FILE
    fi

    if [[ -e $OUTPUT_FILE_RAW ]]; then
	rm $OUTPUT_FILE_RAW
	touch $OUTPUT_FILE_RAW
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
	    HISTNAME=inv_mass_2photon_pi0calib_sector
	    HISTNAME_RAW=inv_mass_2photon_pi0calib_sector_raw

	    root -b -q extract_pi0peak.C\($RUNINDEX,$RUNNUMBER,\"$HISTFILE\",\"$HISTNAME\",$WRITEPLOTS\) | grep "^pi0peak" >> $OUTPUT_FILE
	    root -b -q extract_pi0peak.C\($RUNINDEX,$RUNNUMBER,\"$HISTFILE\",\"$HISTNAME_RAW\",$WRITEPLOTS\) | grep "^pi0peak" >> $OUTPUT_FILE_RAW
	fi

    done

done

exit
