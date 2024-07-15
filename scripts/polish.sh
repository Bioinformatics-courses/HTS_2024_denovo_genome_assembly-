#! /bin/bash

#define directories and file paths
WORKDIR="$HOME/denovo"
ASSEMBLY="$WORKDIR/assembly_output"
SHORT1="$QC_READS/fastp_output/paired_SRR29409522_1.fastq"
SHORT2="$QC_READS/fastp_output/paired_SRR29409522_2.fastq"
POLISH="$WORKDIR/polishing"
MAP="$WORKDIR/mapping"

THREADS=4 #adjust as per your pc caps


mkdir -p $POLISH
mkdir -p $MAP

#First mapping the illumina short reads to the assembly.fasta to generate our BAM file
echo "creating a copy of assembly.fasta in $MAP"
cp -i $ASSEMBLY/assembly.fasta $MAP

echo "now indexing the assembly.fasta file"
bwa index $MAP/assembly.fasta

echo "mapping the illumina short reads to the indexed assembly.fasta. The output is then converted from SAM to BAM and sorted"
bwa mem -t $THREADS $MAP/assembly.fasta $SHORT1 $SHORT2 | samtools view - -Sb | samtools sort - -@ $THREADS -o $MAP/mapped/mapping.bam

echo "indexing mapped output"
bwa index $MAP/mapped/mapping.bam

echo "now polishing the assembled genome using Pilon"
pilon --genome $MAP/assembly.fasta --frags $MAP/mapped/mapping.bam --fix all --output $POLISH/polished

if [ $? -ne 0 ]; then
    echo "polishing process failed"
    exit 1
fi
echo "Post-assembly polishing completed. Output stored in $POLISH"





