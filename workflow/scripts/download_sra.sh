#!/bin/bash

#SBATCH --job-name=download_sra
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=jandr
#SBATCH --time=1:00:00

# Read in args
sra=$1
outdir=$2

module load anaconda; source activate snakemake

# Run fastq-dump
parallel-fastq-dump --sra-id ${sra} --outdir ${outdir} --split-files --gzip 