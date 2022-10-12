#!/bin/bash
# VCF to fasta updated. 

# Arguments
ref=$1
sample_name=$2 # for proper renaming of 'Chromosome'
unfilt_vcf=$3
bed=$4
unmasked_fasta=$5
depth=$6
qual=$7

# Module load tabix

# Name the masked fasta based on the unmasked name.
ppe_masked_fasta=${unmasked_fasta/.fa/_PPEmask.fa}
qfilt_ppe_masked_fasta=${unmasked_fasta/.fa/_qPPEmask.fa}

# Get sample name for correct genotype
samp=$(bcftools query -l ${unfilt_vcf} )

# Index vcf
tabix -p vcf $unfilt_vcf

echo 'creating fasta files for ' $samp

# Check that bcftools is correct version.
bcftools_version=$(bcftools version | head -n1 | grep -Eo '[0-9.]+')
  echo 'bcftools version:' "$bcftools_version"

if [  "$bcftools_version" == 1.9 ] ; then 
  echo 'Incorrect bcftools version'
  exit
fi

# output 1: consensus with no masking (exclude indels).
bcftools consensus --include 'TYPE!="indel"' --fasta-ref ${ref} \
  --sample ${samp} --absent 'N' --missing 'N' ${unfilt_vcf}  | \
  sed "s/>NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome/>$sample_name/g" > ${unmasked_fasta}
    
# seqtk adds a 1 to the sample name, don't use this for renaming. 

# output 2: consensus with ppe masking. 
bcftools consensus --include 'TYPE!="indel"' --mask ${bed} --fasta-ref ${ref} \
  --sample ${samp} --absent 'N' --missing 'N' ${unfilt_vcf}  | \
  sed "s/>NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome/>$sample_name/g" > ${ppe_masked_fasta}
    
echo 'masking with' $bed

# output 3: consensus with ppe masking and quality filters applied.
# Include non-indels with depth >= threshold. For variant sites only (GT==1), filter on QUAL. (No QUAL score for invariant sites.)
bcftools consensus --mask ${bed} --fasta-ref ${ref} \
  --sample ${samp} --include "(TYPE!='indel' & INFO/DP >= $depth) & (QUAL >= $qual | GT == '0')" --absent 'N' --missing 'N' ${unfilt_vcf}  | \
  sed "s/>NC_000962.3 Mycobacterium tuberculosis H37Rv, complete genome/>$sample_name/g" > ${qfilt_ppe_masked_fasta}
    
# Check that masking was correct
filtered_sites=$(grep -v '>' ${qfilt_ppe_masked_fasta} | grep -o 'N' | wc -l )
echo 'applying both the PPE filter and qual/depth filters excluded' $filtered_sites 'sites'