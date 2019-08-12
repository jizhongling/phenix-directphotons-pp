#!/bin/bash
# Function: run failed jobs in pythiaToHisto.sh
# Usage: $0 $(Process)
# Before run this, run in tcsh: echo "\n" >> "${SPIN}/data/pythiaToHisto/aa_badlist.txt"

badlist=(`tail -3 "${SPIN}/data/pythiaToHisto/aa_badlist.txt" | head -2`)
./pythiaToHisto.sh "${badlist[${1}]}"
