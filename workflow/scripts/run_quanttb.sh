#!/bin/bash
# Detect mixed infections in FASTQ files

module load anaconda
source activate quanttb # Need to activate separate python 2.7 environment (https://github.com/AbeelLab/quanttb/issues/7)

p1=$1
p2=$2
output_file=$3

# Detect evidence of mixed infections from FASTQ
quanttb quant -f ${p1} ${p2} -abres -resout -o ${output_file}
