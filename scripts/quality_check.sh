#!/bin/bash

# Directory and parameter definitions
WORKDIR="$HOME/denovo"
RAW_DIR="$WORKDIR/data/raw"
OUT_DIR="$WORKDIR/qc/fastqc_output"
MULTIQC_OUT="$WORKDIR/qc/multiqc_output"
THREADS=4
PAIRED_READ_PREFIX="SRR29409522" #prefix for Illumina read, to be used for Multiqc aggregating --change as needed

# Create output directories if they do not exist
mkdir -p $OUT_DIR
mkdir -p $MULTIQC_OUT

# Loop through each FastQ file and run FastQC
for file in $RAW_DIR/*.{fastq,fastq.gz}
do
    if [[ -e "$file" ]]; then
        echo "Running fastqc on $file"
        fastqc -o $OUT_DIR $file -t $THREADS
    else
        echo "No matching files found in $RAW_DIR"
    fi
done

# Run MultiQC to aggregate FastQC results
echo "Running MultiQC to aggregate fastqc output of paired-end Illumina reads"
multiqc $OUT_DIR/$PAIRED_READ_PREFIX* -o $MULTIQC_OUT

echo "QC check process completed"
