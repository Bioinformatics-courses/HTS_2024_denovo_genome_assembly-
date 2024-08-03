#!/bin/bash
#This scripts runs evaluation on the polished assembly using QUAST

#directory definition
WORKDIR="$HOME/denovo"
POLISH="$WORKDIR/polishing"
EVALUATION_OUT="$WORKDIR/evaluation_output"

#path to files
PILON_POLISHED_FILE="$POLISH/polished.fasta"
REFERENCE_GENOME="$WORKDIR/data/ref/*.fna*"
GFF_FILE="$WORKDIR/data/ref/*.gff*"

THREADS=4 #do adjust as per your pc needs

mkdir -p $EVALUATION_OUT



# Checking if files to evaluate exist
if [ ! -f $GFF_FILE ]; then
    echo "GFF file not found: $GFF_FILE"
    exit 1
fi

if [ ! -f $REFERENCE_GENOME ]; then
    echo "Reference genome file not found: $REFERENCE_GENOME"
    exit 1
fi



# Evaluate the Pilon polished assembly with QUAST
echo -e "\n=============================="
echo "Evaluating the Pilon polished assembly with QUAST"
echo -e "==============================\n"
quast $PILON_POLISHED_FILE -r $REFERENCE_GENOME -g $GFF_FILE -o $EVALUATION_OUT/ --threads $THREADS


echo "Evaluation completed. Results stored in $EVALUATION_OUT"
