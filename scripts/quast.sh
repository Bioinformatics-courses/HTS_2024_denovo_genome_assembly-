#!/bin/bash

#directory definition
WORKDIR="$HOME/denovo"
ASSEMBLY="$WORKDIR/assembly_output"
POLISH="$WORKDIR/polishing"
EVALUATION_OUT="$WORKDIR/evaluation_output"

#path to files
ASSEMBLY_FILE="$ASSEMBLY/assembly.fasta"
PILON_POLISHED_FILE="$POLISH/polished.fasta"
REFERENCE_GENOME="$WORKDIR/data/ref/*.fna"

THREADS=4 #do adjust as per your pc needs

mkdir -p $EVALUATION_OUT

echo -e "\n======================================================="
echo "Both the asasembled genome and the polished genome are going to be evaluated with Quast and stored in seperate output directories"
echo -e "========================================================\n"



# Checking if the assembly and reference genome files exist
if [ ! -f $ASSEMBLY_FILE ]; then
    echo "Assembly file not found: $ASSEMBLY_FILE"
    exit 1
fi

if [ ! -f $REFERENCE_GENOME ]; then
    echo "Reference genome file not found: $REFERENCE_GENOME"
    exit 1
fi


# Evaluate the original assembly with QUAST
echo "Evaluating the original assembly with QUAST"
quast $ASSEMBLY_FILE -r $REFERENCE_GENOME -o $EVALUATION_OUT/quast_before_polishing --threads $THREADS

# Evaluate the Pilon polished assembly with QUAST
echo -e "\n=============================="
echo "Evaluating the Pilon polished assembly with QUAST"
echo -e "==============================\n"
quast $PILON_POLISHED_FILE -r $REFERENCE_GENOME -o $EVALUATION_OUT/quast_after_pilon_polishing --threads $THREADS


echo "Evaluation completed. Results stored in $EVALUATION_OUT"
