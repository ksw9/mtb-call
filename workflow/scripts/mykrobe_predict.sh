#!/bin/bash 
# mykrobe_predict.sh
# run mykrobe predictor. input is a bam and output directory. 

module load anaconda
source activate snakemake

# input
bam=$1
output=$2

# mykrobe predict
mykrobe predict ${bam} tb \
--output ${output}  \
--format csv \
--ploidy haploid \
--seq ${bam} \
--panel 201901	# tb 201901	AMR panel based on first line drugs from NEJM-2018 variants (DOI 10.1056/NEJMoa1800474), and second line drugs from Walker 2015 panel.	NC_000962.3
