#!/bin/bash
# Temporary script to mask PPE regions defined in bed file. 

# Set up environment.
module load seqkit
module load bedtools

# Inputs
unmasked_fa=$1
bed=$2

# Name outputs
renamed_fa=${unmasked_fa/.fa/_rename.fa}
masked_fa=${unmasked_fa/.fa/_masked.fa}

# Get name of fasta. 
sample_name=$(seqkit seq --name $unmasked_fa)
echo $sample_name

# Replace name of fasta for bedtools to work. 
sed "s/>$sample_name/>Chromosome/g" $unmasked_fa > $renamed_fa

# Mask renamed fasta. 
bedtools maskfasta -fi $renamed_fa -bed $bed -fo $masked_fa

# Rename masked fasta in place. 
sed -i "s/>Chromosome/>$sample_name/g" $masked_fa

# Remove temporary file
rm $renamed_fa