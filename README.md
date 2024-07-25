# HTS_2024_Denovo_genome_assembly-

## Introduction
This project demonstrates the denovo genome assembly of the WGS of Lactobacillus plantarum strain JM015 as part of a course requirement. The goal is to assemble the genome using various bioinformatics tools and evaluate the quality of the assembly.

## Setting up our working environment
> **NOTE!**
> 
> This project was demonstrated on a Linux system

Throughout this project, we'll be making use of a couple of tools which will be installed using [conda](https://conda.io/docs/).

Clone the github repo in your home folder and set up the working directory and data directory
```bash
cd
git clone https://github.com/Bioinformatics-courses/HTS_2024_denovo_genome_assembly-.git
mv HTS_2024_denovo_genome_assembly- denovo
cd denovo
mkdir -p data/raw

```

> **CODE BREAKDOWN**
> 
> - **`cd`** - changing to the $HOME directory (if not already there)
> - **`mv HTS_2024_denovo_genome_assembly- denovo`** - renames the cloned dir to "denovo" (this dir is set in the sripts as the main working directory)
> - **`cd denovo`** - takes us to the denovo directory
> - **`mkdir -p data/raw`** - creates data/raw where we will keep our fastq raw files
> 
Next we create our conda environment from the yml file provided in the env folder and activate the environment
```bash
conda env create -f env/denovo-conda.yml
conda activate denovo

```

## The Data
The datasets used for this assemby process are WGS of Lactiplantibacillus plantarum (formerly Lactobacillus plantarum). Since we are going to be doing hybrid assembly, we are using both short and long reads. The datasets here were sequenced from the same biological sample and consits of :
- A long read PACBIO_SMRT run with accession number [SRR29409521](https://www.ncbi.nlm.nih.gov/sra/SRX24922988[accn])
- Paired-end ILLUMINA short reads with accession number [SRR29409522](https://www.ncbi.nlm.nih.gov/sra/SRX24922987[accn])

They can be sourced from the links provided above, but better still, can be obtained using sra-toolkit
```bash
sudo apt-get install sra-toolkit 
prefetch SRR29409521 SRR29409522 -O data/
fasterq-dump data/SRR29409521 --outdir data/raw
fasterq-dump data/SRR29409522 --outdir data/raw

```
We'll end up with 3 files in the data/raw subdirectory of our working directory, SRR29409521.fastq, SRR29409522_1.fastq and SRR29409522_2.fastq, which are our long read and paired-end short reads respectively.

## Quality check
The first script we're going to run is the quality_check.sh. This script runs fastqc on the fastq files in data/raw and generates a fastqc report. It also aggregates the report of the paired-end read with Multiqc.

```bash
cd scripts
sudo chmod +x quality_check.sh
./quality_checks.sh
```

> **CODE BREAKDOWN**
> 
> - **`cd scripts`** - this moves us into our scripts directory
> - **`sudo chmod +x quality_check.sh`** - grants execute permission to the script
> - **`./quality_checks.sh`** - executes the script
>

Taking a look at our script (properly commented of course), the script sets up the working directories and creates our output directory quality_check_reports/ after which it loops over the files in dir/raw and runs fastqc on them. It also aggregates the fastqc report of the paired-end reads with MultiQC.

Visualizing the output of our long read quality check (denovo/quality_check_report/fastqc_output/SRR29409521_fastqc.html). The report shows overall high per base sequence quality, with median scores above 90 and only a low average dip to 66 towards the end.
<center><img src="_static/long_read.png" width="90%"></center>

Next, we examine the MultiQC report of our paired-end read (denovo/quality_check_report/multiqc_output/multiqc_report.html). We also see overall good reports here with the mean quality score generally hovering around the 35 mark.
<center><img src="_static/fastqc_per.png" width="90%"></center>

## Quality Control
Before we do proceed to assembly, we have to do some quality control. To do this, we use the script 






