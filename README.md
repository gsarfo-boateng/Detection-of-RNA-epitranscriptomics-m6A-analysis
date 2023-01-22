# Detection-of-RNA-epitranscriptomics-m6A-analysis
RNA modification analysis python3, Bash and Snakemake
George Boateng

the python script are for the analysis of m6ACE that were processed by Dr. Sho Goh
Bash scripts are for easy to proccess fastq files for metanalysis.
Snakemake is a Bioinformatics tool for managing a workflow. This tool proves valuable when analyzing a large amount of data with multiple tools. This script was made as a learning tool for workflow manager. There is also Nextflow to manage large analysis workflows. Here, Snakemake was used to run everything that is usually run on Linux with RNA-Seq Analyses (here is my long winded version of an RNA-seq analysis on Mice p53 gene mutation).


The genome and the gtf files were downloaded and an index was created of the genome using the dm6.fasta and gtf files from NCBI.

Conda Tools Used:

FastQC
fastp
STAR
featureCounts (conda Subread)
Bowtie2
bwa
Samtools
MultiQC
Running

Python scripts for m6ACE have defined directories. Move samples that you wish to analysis into that directories.
Run snakemake in the same folder as the snakefile used snakemake -j 80 which tells Snakemake to use 80 cores. Snakemake if not given a file name searches current directory for a file named snakefile.
