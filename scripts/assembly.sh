#!/bin/bash

# Script to perform hybrid assembly using Unicycler

# Set directories
WORKDIR="$HOME/denovo"
QC_READS="$WORKDIR/data/QC"
ASSEMBLY_OUT="$WORKDIR/assembly_output"

#paths to the processed reads
SHORT1="$QC_READS/fastp_output/paired_SRR29409522_1.fastq"
SHORT2="$QC_READS/fastp_output/paired_SRR29409522_2.fastq"
UNPAIRED_SHORT="$QC_READS/fastp_output/unpaired_SRR29409522_1.fastq"
LONG="$QC_READS/filtlong_output/filtered_long_read.fastq"

THREADS=4 #default unicycler usage, adjust as necessary for faster assembly
# Creating the output directory if it doesn't exist
mkdir -p $ASSEMBLY_OUT


# Run Unicycler for hybrid assembly
echo "Running Unicycler for hybrid assembly"
unicycler -1 $SHORT1 -2 $SHORT2 -s $UNPAIRED_SHORT -l $LONG -o $ASSEMBLY_OUT --keep 2 --thread $THREADS

if [ $? -ne 0 ]; then
    echo "Error: Assembly failed"
    exit 1
fi
echo "Hybrid assembly completed. Output stored in $ASSEMBLY_OUT"

