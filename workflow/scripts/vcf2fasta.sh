#!/bin/bash
# VCF to fasta

# read from command line.
ref=$1
vcf=$2
sample_name=$3
bed=$4
fasta=$5

# Get sample name for correct genotype
samp=$(bcftools query -l ${vcf} )

# Index vcf
tabix -p -f vcf $vcf

# Make consensus masked sequence & rename fasta header. 
bcftools consensus --include 'TYPE!="indel"' --mask ${bed} --fasta-ref ${ref} \
  --sample ${samp} --absent 'N' --missing 'N' ${vcf}  | \
  sed "s/>NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome/>$sample_name/g" > ${fasta}
    
# seqtk adds a 1 to the sample name, don't use.