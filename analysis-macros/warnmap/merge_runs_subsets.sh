#! /bin/zsh

DATASETS=(
    Run13pp510MinBias
)

INPUT_DIR_BASE="/phenix/spin3/nfeege/taxi/"
OUTPUT_DIR="./warnmap-data/"

for DATASET in $DATASETS; do

    echo "Processing data set $DATASET ..."

    INPUT_DIR="$INPUT_DIR_BASE/$DATASET"

    [[ ! -d $INPUT_DIR ]] && echo "ERROR: Directory $INPUT_DIR does not exist. EXIT." && exit

    TAXI_DATA_DIR=`awk '$2=="Run_DirectPhotonPP_Warnmap" {print $1}' $INPUT_DIR/taxilist.txt`

    INPUT_DIR="$INPUT_DIR/$TAXI_DATA_DIR/data/"

    echo "Merging files from input directory $INPUT_DIR to $OUTPUT_FILE ... "

    OUTPUT_LIST_SET1="${OUTPUT_DIR}/WarnmapData_${DATASET}_RunRange1_List.txt"
    OUTPUT_LIST_SET2="${OUTPUT_DIR}/WarnmapData_${DATASET}_RunRange2_List.txt"
    OUTPUT_LIST_SET3="${OUTPUT_DIR}/WarnmapData_${DATASET}_RunRange3_List.txt"
    OUTPUT_LIST_SET4="${OUTPUT_DIR}/WarnmapData_${DATASET}_RunRange4_List.txt"
    OUTPUT_LIST_SET5="${OUTPUT_DIR}/WarnmapData_${DATASET}_RunRange5_List.txt"

    ls $INPUT_DIR/WarnmapData*.root | head -200  | tail -200 > $OUTPUT_LIST_SET1
    ls $INPUT_DIR/WarnmapData*.root | head -400  | tail -200 > $OUTPUT_LIST_SET2
    ls $INPUT_DIR/WarnmapData*.root | head -600  | tail -200 > $OUTPUT_LIST_SET3
    ls $INPUT_DIR/WarnmapData*.root | head -800  | tail -200 > $OUTPUT_LIST_SET4
    ls $INPUT_DIR/WarnmapData*.root | head -1000  | tail -200 > $OUTPUT_LIST_SET5

    FILELIST_SET1=("${(@f)$(ls $INPUT_DIR/WarnmapData*.root | head -200  | tail -200)}")
    FILELIST_SET2=("${(@f)$(ls $INPUT_DIR/WarnmapData*.root | head -400  | tail -200)}")
    FILELIST_SET3=("${(@f)$(ls $INPUT_DIR/WarnmapData*.root | head -600  | tail -200)}")
    FILELIST_SET4=("${(@f)$(ls $INPUT_DIR/WarnmapData*.root | head -800  | tail -200)}")
    FILELIST_SET5=("${(@f)$(ls $INPUT_DIR/WarnmapData*.root | head -1000 | tail -200)}")

    OUTPUT_FILE_SET1="${OUTPUT_DIR}/WarnmapData_${DATASET}_RunRange1.root"
    OUTPUT_FILE_SET2="${OUTPUT_DIR}/WarnmapData_${DATASET}_RunRange2.root"
    OUTPUT_FILE_SET3="${OUTPUT_DIR}/WarnmapData_${DATASET}_RunRange3.root"
    OUTPUT_FILE_SET4="${OUTPUT_DIR}/WarnmapData_${DATASET}_RunRange4.root"
    OUTPUT_FILE_SET5="${OUTPUT_DIR}/WarnmapData_${DATASET}_RunRange5.root"

    hadd $OUTPUT_FILE_SET1 $FILELIST_SET1
    hadd $OUTPUT_FILE_SET2 $FILELIST_SET2
    hadd $OUTPUT_FILE_SET3 $FILELIST_SET3
    hadd $OUTPUT_FILE_SET4 $FILELIST_SET4
    hadd $OUTPUT_FILE_SET5 $FILELIST_SET5

    echo "DONE."

done

exit
