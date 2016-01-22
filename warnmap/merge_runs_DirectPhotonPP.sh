#! /bin/zsh

DATASETS=(
    Run13pp510MinBias
    Run13pp510ERT
)

INPUT_DIR_BASE="/phenix/spin3/nfeege/taxi/"
OUTPUT_DIR="./warnmap-data/"

for DATASET in $DATASETS; do

    echo "Processing data set $DATASET ..."

    OUTPUT_FILE="${OUTPUT_DIR}/DirectPhotonPP_${DATASET}.root"

    INPUT_DIR="$INPUT_DIR_BASE/$DATASET"

    [[ ! -d $INPUT_DIR ]] && echo "ERROR: Directory $INPUT_DIR does not exist. EXIT." && exit

    TAXI_DATA_DIR=`awk '$2=="Run_DirectPhotonPP" {print $1}' $INPUT_DIR/taxilist.txt`

    INPUT_DIR="$INPUT_DIR/$TAXI_DATA_DIR/data/"

    echo "Merging files from input directory $INPUT_DIR to $OUTPUT_FILE ... "

    ## LIMIT hadd to 1000 files at a time because of "Too many open files" error
    NFILES=$( ls $INPUT_DIR/DirectPhotonPP*.root | wc -l )
    echo "Number of files to merge: $NFILES"

    if [[ $NFILES -le 1000 ]]; then

	FILELIST=("${(@f)$(ls $INPUT_DIR/DirectPhotonPP*.root)}")

#        haddPhenix $OUTPUT_FILE $FILELIST
	hadd $OUTPUT_FILE $FILELIST

    elif [[ $NFILES -le 2000 ]]; then

	echo "File list between 1000 and 2000 entries long, need to split..."

	N_SUB1=1000;
	let "N_SUB2 = $NFILES - $N_SUB1"

	echo $NFILES
	echo $N_SUB1
	echo $N_SUB2

	FILELIST_SUB1=("${(@f)$(ls $INPUT_DIR/DirectPhotonPP*.root | head -${N_SUB1})}")
	FILELIST_SUB2=("${(@f)$(ls $INPUT_DIR/DirectPhotonPP*.root | tail -${N_SUB2})}")

	OUTPUT_FILE_SUB1="${OUTPUT_DIR}/DirectPhotonPP_${DATASET}_Temporary_Sub1.root"
	OUTPUT_FILE_SUB2="${OUTPUT_DIR}/DirectPhotonPP_${DATASET}_Temporary_Sub2.root"

	hadd $OUTPUT_FILE_SUB1 $FILELIST_SUB1
	hadd $OUTPUT_FILE_SUB2 $FILELIST_SUB2

	hadd $OUTPUT_FILE $OUTPUT_FILE_SUB1 $OUTPUT_FILE_SUB2

	rm $OUTPUT_FILE_SUB1 $OUTPUT_FILE_SUB2

    else

	echo "File list longer than 2000 entries, stop processing here."
	return

    fi

    echo "DONE."

done

exit
