#! /bin/bash

#define directories and file paths
WORKDIR="$HOME/denovo"
ASSEMBLY="$WORKDIR/assembly_output"
SHORT1="$WORKDIR/data/QC/fastp_output/paired_SRR29409522_1.fastq"
SHORT2="$WORKDIR/data/QC/fastp_output/paired_SRR29409522_2.fastq"
POLISH="$WORKDIR/polishing"
MAP="$WORKDIR/mapping"

THREADS=4 #adjust as per your pc caps


mkdir -p $POLISH
mkdir -p $MAP
mkdir -p $MAP/mapped

#First mapping the illumina short reads to the assembly.fasta to generate our BAM file
echo -e "\n=============================="
echo "creating a copy of assembly.fasta in $MAP"
echo -e "==============================\n"

cp -i $ASSEMBLY/assembly.fasta $MAP

echo -e "\n=============================="
echo "now indexing the assembly.fasta file"
echo -e "==============================\n"

bwa index $MAP/assembly.fasta

echo -e "\n=============================="
echo "mapping the illumina short reads to the indexed assembly.fasta. The output is then converted from SAM to BAM and sorted"
echo -e "==============================\n"

bwa mem -t $THREADS $MAP/assembly.fasta $SHORT1 $SHORT2 | samtools view - -Sb | samtools sort - -@ $THREADS -o $MAP/mapped/mapping.bam

echo -e "\n=============================="
echo "indexing sorted BAM output with samtools"
echo -e "==============================\n"
samtools index $MAP/mapped/mapping.bam -@ $THREADS

echo -e "\n=============================="
echo "now polishing the assembled genome using Pilon"
echo -e "==============================\n"
pilon --genome $MAP/assembly.fasta --frags $MAP/mapped/mapping.bam --fix all --changes --output $POLISH/polished

if [ $? -ne 0 ]; then
    echo "polishing process failed"
    exit 1
fi

echo -e "\n=============================="
echo "Post-assembly polishing completed. Output stored in $POLISH"
echo -e "==============================\n"





