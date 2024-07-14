#!/bin/bash

# Directory and parameter definitions
WORKDIR="$HOME/denovo"
RAW_DIR="$WORKDIR/data/raw"
OUT_DIR="$WORKDIR/quality_check_reports/fastqc_output"
MULTIQC_OUT="$WORKDIR/quality_check_reports/multiqc_output"

THREADS=4 #adjust as needed
PAIRED_READ_PREFIX="SRR29409522" #prefix for Illumina read, to be used for Multiqc aggregating --change as needed

# Create output directories if they do not exist
mkdir -p $OUT_DIR
mkdir -p $MULTIQC_OUT

# Looping through each fastq || fastq.gz file and run fastqc
for file in $RAW_DIR/*.{fastq,fastq.gz}
do
    if [[ -e "$file" ]]; then
        echo "Running fastqc on $file"
        fastqc -o $OUT_DIR $file -t $THREADS
    else
        echo "No matching files found in $RAW_DIR"
    fi
done

# Run MultiQC to aggregate FastQC results for short pair-end reads
echo "Running MultiQC to aggregate fastqc output of paired-end Illumina reads"
multiqc $OUT_DIR/$PAIRED_READ_PREFIX* -o $MULTIQC_OUT

echo "Quality check process completed"
echo "Output files can be found in "$WORKDIR/quality_chech_reports"
