#!/bin/zsh

# OBSOLETE - Use PYTHON script instead!!!

FILELIST=(
    'warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange0.txt'
    'warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange1.txt'
    'warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange2.txt'
    'warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange3.txt'
    'warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange4.txt'
)

# summary arranged by warnmap
for FILE in $FILELIST; do
    echo ""
    echo "BEGIN summary file $FILE"

    for SECTOR in {0..7}
    do
	echo "Hot towers Sector $SECTOR: $(grep ^${SECTOR} $FILE | grep 50$ | wc -l)"
    done

    echo "END summary file $FILE"
done

# summary arranged by sector
for SECTOR in {0..7}
do
    echo ""
    echo "BEGIN summary sector $SECTOR"

    for FILE in $FILELIST; do
	echo "Hot towers file $(echo $FILE | cut -d_ -f5 | cut -d. -f1): $(grep ^${SECTOR} $FILE | grep 50$ | wc -l)"
    done

    echo "END summary sector $SECTOR"
done

