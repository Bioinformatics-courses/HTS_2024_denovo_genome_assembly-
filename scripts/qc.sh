#!/bin/bash/

#directory and parameter definitions
WORKDIR="/home/denovo"
RAW_DIR="$WORKDIR/raw"
OUT_DIR="$WORKDIR/qc/fastqc_output"
THREADS=4
#including multiqc aggregation
MULTIQC_OUT="$WORKDIR/qc/multiqc_output"

mkdir -p $OUT_DIR

for file in $RAW_DIR/*.{fastq,fastq.gz,fq,fq.gz}
do
	if [ -e "$file"]; then
		echo "Running fastqc on $file"
		fastqc -o $OUT_DIR $file -t $THREADS
	else
		echo "No matching files found in $RAW_DIR"
	fi
done


