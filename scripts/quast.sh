#!/bin/bash

#directory definition
WORKDIR="$HOME/denovo"
ASSEMBLY="$WORKDIR/assembly_output"
POLISH="$WORKDIR/polishing/polished"

THREADS=4 #do adjust as per your pc needs

echo -e "\n======================================================="
echo "Both the asasembled genome and the polished genome are going to be evaluated with Quast and stored in seperate directories"
echo -e "========================================================\n"



