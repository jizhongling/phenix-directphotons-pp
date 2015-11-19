#!/bin/zsh

# OBSOLETE - use PYTHON script instead!!!

# flow: check all input files for all sector/y/z combination. If any status if 50, use 50. If any status 40, use 40. If any status 20, use 20. If any status 10, use 10.
# Otherwise, use 0.

#Loop over lines in first file and check status in other files.
#warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange0.txt
#warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange1.txt
#warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange2.txt
#warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange3.txt
#warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange4.txt

grep 50$ warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_erange*.txt | cut -d: -f2 | sort -u > warnmap_output/Warnmap_Run13pp510MinBias_mergeruns_mergeerange.txt
