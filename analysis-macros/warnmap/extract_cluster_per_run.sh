#! /bin/zsh

DATASETS=(
#    Run13pp510MinBias
    Run13pp510ERT
)

INPUT_DIR_BASE="/phenix/spin3/nfeege/taxi/"
OUTPUT_DIR="./warnmap_data/"

for DATASET in $DATASETS; do

    OUTPUT_FILE="allruns_${DATASET}_cluster_temp.txt"

    if [[ -e $OUTPUT_FILE ]]
    then
	rm $OUTPUT_FILE
    fi

    touch $OUTPUT_FILE

    echo "Processing data set $DATASET ..."

    INPUT_DIR="$INPUT_DIR_BASE/$DATASET"

    [[ ! -d $INPUT_DIR ]] && echo "ERROR: Directory $INPUT_DIR does not exist. EXIT." && exit

    TAXI_DATA_DIR=`awk '$2=="Run_DirectPhotonPP_Warnmap" {print $1}' $INPUT_DIR/taxilist.txt`

    INPUT_DIR="$INPUT_DIR/$TAXI_DATA_DIR/data/"

    echo "Checking number of cluster in files from input directory $INPUT_DIR ... "

    FILELIST=`ls $INPUT_DIR | grep "WarnmapData"`

    for FILE in $INPUT_DIR/WarnmapData*.root; do

	#echo "Check $FILE"

	I_RUN=$( echo $FILE | cut -d \- -f2 | cut -d. -f1 )
	I_CLUS=$( root -b -q extract_cluster_per_run.C\(\"$FILE\"\) | grep "Cluster / event in run:" )

	echo "$I_RUN $I_CLUS" >> $OUTPUT_FILE

    done

    echo "DONE."

done

exit
