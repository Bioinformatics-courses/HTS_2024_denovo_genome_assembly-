#!/bin/bash

# Directory and parameter definitions
WORKDIR="$HOME/denovo"
RAW_DIR="$WORKDIR/data/raw"
FILTLONG_OUT="$WORKDIR/data/QC/filtlong_output"
FASTP_OUT="$WORKDIR/data/QC/fastp_output"
REPORT="$FASTP_OUT/report"

#direct path to reads
LONG_READ="$RAW_DIR/SRR29409521.fastq"
SHORT_READ1="$RAW_DIR/SRR29409522_1.fastq"
SHORT_READ2="$RAW_DIR/SRR29409522_2.fastq"

THREADS=4  #adjust as needed

#create output directories if they don't exist
mkdir -p $FILTLONG_OUT
mkdir -p $FASTP_OUT
mkdir -p $REPORT

echo "Running Filtlong on $LONG_READ"
filtlong --min_length 1000 --keep_percent 95 $LONG_READ > $FILTLONG_OUT/filtered_long_read.fastq
if [ $? -ne 0 ]; then
    echo "Error: Filtlong failed to run on $LONG_READ"
    exit 1
fi

echo "Running fastp on $SHORT_READ1 and $SHORT_READ2"
fastp --in1 $SHORT_READ1 \
      --in2 $SHORT_READ2 \
      --out1 $FASTP_OUT/paired_SRR29409522_1.fastq \
      --out2 $FASTP_OUT/paired_SRR29409522_2.fastq \
      --unpaired1 $FASTP_OUT/unpaired_SRR29409522_1.fastq \
      --unpaired2 $FASTP_OUT/unpaired_SRR29409522_2.fastq \
      --thread $THREADS \
      --html $REPORT/fastp.html \
      --json $REPORT/fastp.json
if [ $? -ne 0 ]; then
    echo "Error: fastp failed to run on $SHORT_READ1 and $SHORT_READ2"
    exit 1
fi

echo "Quality control process completed. Cleaned Filtlong output stored in $FILTLONG_OUT and cleaned fastp output in $FASTP_OUT"
