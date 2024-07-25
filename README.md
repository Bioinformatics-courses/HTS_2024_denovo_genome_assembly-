# HTS_2024_Denovo_genome_assembly-

## Introduction
This project demonstrates the denovo genome assembly of the WGS of Lactobacillus plantarum strain JM015 as part of a course requirement. The goal is to assemble the genome using various bioinformatics tools and evaluate the quality of the assembly. We're going to be 

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
Before we do proceed to assembly, we have to do some quality control. Our datasets from our quality check do not look to be in need of trimming, or adapter removal. So instead, we're just going to run some basic quality control on our short paired-end reads using the default qc parameters of the tool [fastp](https://github.com/OpenGene/fastp) in our script. It’s easy to use and lives up to its name (very fast).

For our PacBio long read data, we're also going to be running very light QC, using the tool [Filtlong](https://github.com/rrwick/Filtlong) in our script. Filtlong is set in the script to remove any reads shorter than 1 kbp and also exclude the worst 5% of reads. The reason for this light QC is that the data is of high quality, so we just want to filter out the worst of the read. Filtlong prefers longer reads and aggressively filters out shorter ones. This generally benefits genome assembly but can be detrimental for small plasmids, hence why we only choose to exclude reads shorter than 1kpb.

To run the quality_control.sh script, we navigate to our scripts folder, give execute permission to quality_control.sh and run the file.

```bash
cd scripts # assuming we're in our working folder
sudo chmod +x quality_control.sh
./quality_control.sh
```
This will generate 4 outputs for the short reads: 2 paired and 2 unpaired reads for our fastp QC, and a filtered long read file from the Filtlong QC, all located in the corresponding subdirectories under data/QC.
The unpaired reads from the fastpc QC are very small compared to the the paired reads and that is very good. Taking a look at the fastp report generated, we can see that not much changed.

## Assembly
For the assembly, we'll be running the assembly.sh script in the scripts folder. The primary tool in this script is [Unicycler](https://github.com/rrwick/Unicycler). Unicycler is an assembly pipeline for bacterial genomes and comes highly recommended for hybrid assembly of bacteria.
Unicycler uses [SPAdes] to produce graphs, which are made by performing a de Bruijn graph assembly with a range of different k-mer sizes. After generating a short-read assembly graph, Unicycler then uses the long reads to scaffold the graph to completion. It performs best when the short reads are deep and have even coverage.
Unicycler hybrid assembly is run by passing our paired-end short reads, and our unpaired(if available) and the final secret sauce, the long read, to Unicycler, which is exactly what we have in out script. We run the script using:

```bash
cd scripts # assuming we're in our working folder
sudo chmod +x assembly.sh
./assembly.sh
```

Unicycler takes a long time to run hybrid assembly; with 4 threads given on my i5 5th gen cpu with 8gbs of ram, it took 11 hours to complete the hybrid assembly! The output is stored in assembly_output/ where we have the most important files: assembly.fasta (our assembled genome), assembly.gfa (our assembly graph) and unicycler.log (log file).

## Polishing




