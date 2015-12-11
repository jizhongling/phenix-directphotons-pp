#! /bin/zsh

DATASETS=(
#    Run9pp200MinBias
#    Run9pp200ERT
#    Run9pp500MinBias
#    Run9pp500ERT
    Run13pp510MinBias
    Run13pp510ERT
)

INPUT_DIR_BASE="/phenix/spin3/nfeege/taxi/"
OUTPUT_DIR="./warnmap-data/"

for DATASET in $DATASETS; do

    echo "Processing data set $DATASET ..."

    OUTPUT_FILE="${OUTPUT_DIR}/WarnmapData_${DATASET}.root"

    INPUT_DIR="$INPUT_DIR_BASE/$DATASET"

    [[ ! -d $INPUT_DIR ]] && echo "ERROR: Directory $INPUT_DIR does not exist. EXIT." && exit

    TAXI_DATA_DIR=`awk '$2=="Run_DirectPhotonPP_Warnmap" {print $1}' $INPUT_DIR/taxilist.txt`

    INPUT_DIR="$INPUT_DIR/$TAXI_DATA_DIR/data/"

    echo "Merging files from input directory $INPUT_DIR to $OUTPUT_FILE ... "

    ## LIMIT hadd to first 1000 files in folder because of "Too many open files" error
    FILELIST=("${(@f)$(ls $INPUT_DIR/WarnmapData*.root | head -1000)}")

#    haddPhenix $OUTPUT_FILE $FILELIST
    hadd $OUTPUT_FILE $FILELIST

    echo "DONE."

done

exit
