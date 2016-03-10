#! /bin/zsh

ALLOW_SKIP_RUNS=false

DATASETS=(
    Run13pp510MinBias
    Run13pp510ERT
)

INPUT_DIR_BASE="/phenix/spin3/nfeege/taxi/"
OUTPUT_DIR="./data/"

RUNLISTFILE="../runqa/Run13pp510_RunQuality.txt"

IFS=$'\n' RUNLIST=($(cat $RUNLISTFILE | grep 0$ | cut -d\  -f1))
NRUNS=$(cat $RUNLISTFILE | grep 0$ | wc -l)

echo "Using Run List $RUNLISTFILE"

TEMPFILE=templist.txt

for DATASET in $DATASETS; do

    echo "Processing data set $DATASET ..."

    if [[ -e $TEMPFILE ]]; then
	rm $TEMPFILE
	touch $TEMPFILE
    fi

    OUTPUT_FILE="${OUTPUT_DIR}/DirectPhotonPP-${DATASET}.root"

    INPUT_DIR="$INPUT_DIR_BASE/$DATASET"

    [[ ! -d $INPUT_DIR ]] && echo "ERROR: Directory $INPUT_DIR does not exist. EXIT." && exit

    TAXI_DATA_DIR=`awk '$2=="Run_DirectPhotonPP" {print $1}' $INPUT_DIR/taxilist.txt`

    INPUT_DIR="$INPUT_DIR/$TAXI_DATA_DIR/data/"

    echo "Merging files from input directory $INPUT_DIR to $OUTPUT_FILE ... "

    ## LIMIT hadd to 1000 files at a time because of "Too many open files" error
    if [[ $NRUNS -le 1000 ]]; then

	for RUN in $RUNLIST; do
#	    if [[ $ALLOW_SKIP_RUNS="true" ]]; then
#		if [[ -e $INPUT_DIR/DirectPhotonPP-${RUN}.root ]]; then
#		    echo "$INPUT_DIR/DirectPhotonPP-${RUN}.root" >> $TEMPFILE
#		fi
#	    else
		echo "$INPUT_DIR/DirectPhotonPP-${RUN}.root" >> $TEMPFILE
#	    fi
	done

	FILELIST=("${(@f)$(cat $TEMPFILE)}")
	echo $FILELIST
	hadd $OUTPUT_FILE $FILELIST
    else
	echo "More than 1000 runs, stop processing. ERROR."
	exit
    fi

    echo "DONE."

done

exit
