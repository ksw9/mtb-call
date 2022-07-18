#!/bin/bash
# Variant quality score recalibration (VQSR) with GATK. 
# Takes (1) reference genome and (2) VCF file and (3) Output VQSR VCF file. 
# By default, this uses the variants with QUAL > mean QUAL as "truth set" for training. 

module load anaconda
source activate snakemake

# read from command line
ref=$1
vcf=$2
vqsr_vcf=$3
train_vcf=$4

# get basename
base=${vcf%.vcf*}

# Set tranche filter level: this defines the sensitivity for the "truth variants."
ts_filter=99.0

# If training set doesn't exist, merge VCFs and create one.
if [[ -z ${train_vcf} ]]; then
  echo 'no training set provided, creating internal training set'
  
  train_vcf=${vcf/.vcf.gz/_train.vcf.gz}
  
  # get mean QUAL, excluding QUAL of invariant sites, to select variants for training
  qual=$(bcftools filter -i 'TYPE == "SNP"'  ${vcf}  | bcftools query  -f '%QUAL\n' | grep -v inf | \
  awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' )

  echo $qual
  echo 'mean qual of variant sites': $qual

  # Create truth set by selecting only high qual variants (QUAL > mean) among variant sites.
  bcftools view --types snps ${vcf} | bcftools filter -e 'QUAL == inf' | bcftools filter -i " QUAL > $qual " -O z -o ${train_vcf}

  # Index vcf.
  tabix -f -p vcf ${train_vcf}
fi

# Recalibrate variants. 
 
if gatk VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource:truePos,known=false,training=true,truth=true,prior=15 ${train_vcf} \
-an DP \
-an QD \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-an SOR \
-an MQ \
-mode SNP \
--max-gaussians 2 \
--output ${base}.recal \
--tranches-file ${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0 

then echo 'VQSR succeeded with MQRankSum'

   # Apply VQSR.
    gatk ApplyVQSR \
    -R ${ref} \
    -mode SNP \
    --variant ${vcf}  \
    --recal-file ${base}.recal  \
    --tranches-file ${base}".tranches"  \
    --truth-sensitivity-filter-level ${ts_filter} \
    --output ${vqsr_vcf/.vcf.gz/.vcf}

  # Test if this works in calculating VQSLOD. If the VQSLOD score is Nan/inf, rerun, setting lowering number of Gaussians.
  echo 'VQSLOD contains inf/Nan at: ' $(bcftools query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

  # If VQSR does not run, remove MQRankSum.
  failedSites=$(bcftools query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

  if [ "$failedSites" -eq "0" ]; then
    echo 'Apply variants succeeded with MQRankSum'
  fi

# If failed, remove MQRankSum which may not have sufficient variation.
elif gatk VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource:truePos,known=false,training=true,truth=true,prior=15 ${train_vcf} \
-an DP \
-an QD \
-an ReadPosRankSum \
-an FS \
-an SOR \
-an MQ \
-mode SNP \
--max-gaussians 2 \
--output ${base}.recal \
--tranches-file ${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0  

then echo 'VQSR w/o MQRankSum'

   # Apply VQSR.
    gatk ApplyVQSR \
    -R ${ref} \
    -mode SNP \
    --variant ${vcf}  \
    --recal-file ${base}.recal  \
    --tranches-file ${base}".tranches"  \
    --truth-sensitivity-filter-level ${ts_filter} \
    --output ${vqsr_vcf/.vcf.gz/.vcf}

  # Test if this works in calculating VQSLOD. If the VQSLOD score is Nan/inf, rerun, setting lowering number of Gaussians.
  echo 'VQSLOD contains inf/Nan at: ' $(bcftools query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

  # If VQSR does not run, remove MQRankSum.
  failedSites=$(bcftools query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

  if [ "$failedSites" -eq "0" ]; then
    echo 'Apply variants succeeded removing MQRankSum'
  
  
# If failed, remove ReadPosRankSum which may not have sufficient variation.
elif gatk VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource:truePos,known=false,training=true,truth=true,prior=15 ${train_vcf} \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-mode SNP \
--max-gaussians 2 \
--output ${base}.recal \
--tranches-file ${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0  

then echo 'VQSR w/o MQRankSum & ReadPosRankSum'

   # Apply VQSR.
    gatk ApplyVQSR \
    -R ${ref} \
    -mode SNP \
    --variant ${vcf}  \
    --recal-file ${base}.recal  \
    --tranches-file ${base}".tranches"  \
    --truth-sensitivity-filter-level ${ts_filter} \
    --output ${vqsr_vcf/.vcf.gz/.vcf}

  # Test if this works in calculating VQSLOD. If the VQSLOD score is Nan/inf, rerun, setting lowering number of Gaussians.
  echo 'VQSLOD contains inf/Nan at: ' $(bcftools query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

  # If VQSR does not run, remove MQRankSum.
  failedSites=$(bcftools query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

  if [ "$failedSites" -eq "0" ]; then
    echo 'Apply variants succeeded removing MQRankSum & ReadPosRankSum'
  fi
  
# If failed, maxGaussians=1.
elif gatk VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource:truePos,known=false,training=true,truth=true,prior=15 ${train_vcf} \
-an DP \
-an QD \
-an FS \
-an SOR \
-an MQ \
-mode SNP \
--max-gaussians 1 \
--output ${base}.recal \
--tranches-file ${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0  

  then echo 'VQSR w/ 1 Gaussian '

   # Apply VQSR.
    gatk ApplyVQSR \
    -R ${ref} \
    -mode SNP \
    --variant ${vcf}  \
    --recal-file ${base}.recal  \
    --tranches-file ${base}".tranches"  \
    --truth-sensitivity-filter-level ${ts_filter} \
    --output ${vqsr_vcf/.vcf.gz/.vcf}

  # Test if this works in calculating VQSLOD. If the VQSLOD score is Nan/inf, rerun, setting lowering number of Gaussians.
  echo 'VQSLOD contains inf/Nan at: ' $(bcftools query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

  # If VQSR does not run, remove MQRankSum.
  failedSites=$(bcftools query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

  if [ "$failedSites" -eq "0" ]; then
    echo 'Apply variants succeeded with maxGaussian=1'
  fi

# If failed, remove MQ, reset maxGaussians to 2.
elif gatk VariantRecalibrator \
-R ${ref} \
--variant ${vcf} \
--resource:truePos,known=false,training=true,truth=true,prior=15 ${train_vcf} \
-an DP \
-an QD \
-an FS \
-an SOR \
-mode SNP \
--max-gaussians 2 \
--output ${base}.recal \
--tranches-file ${base}".tranches"  \
-tranche 100.0 -tranche 99.0 -tranche 96.0 -tranche 93.0 -tranche 90.0  

  then echo 'VQSR w/ 1 Gaussian '

   # Apply VQSR.
    gatk ApplyVQSR \
    -R ${ref} \
    -mode SNP \
    --variant ${vcf}  \
    --recal-file ${base}.recal  \
    --tranches-file ${base}".tranches"  \
    --truth-sensitivity-filter-level ${ts_filter} \
    --output ${vqsr_vcf/.vcf.gz/.vcf}

  # Test if this works in calculating VQSLOD. If the VQSLOD score is Nan/inf, rerun, setting lowering number of Gaussians.
  echo 'VQSLOD contains inf/Nan at: ' $(bcftools query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

  # If VQSR does not run, remove MQRankSum.
  failedSites=$(bcftools query  -f '%INFO/VQSLOD\n' ${base}_vqsr.vcf | grep '\inf$\|\Nan$|\.$' | wc -l)

  if [ "$failedSites" -eq "0" ]; then
    echo 'Apply variants succeeded with no MQ'
  fi


fi
fi
# Tabix index and zip all VCF files for vcf-merge to work
bgzip -f -c ${vqsr_vcf/.vcf.gz/.vcf} > ${vqsr_vcf}
tabix -f -p vcf ${vqsr_vcf}

# Remove intermediate files.
rm ${base}.recal
rm ${base}.recal.idx
#rm ${train_vcf}
#rm ${train_vcf}.tbi
rm ${base}".tranches"
#rm ${base}_vqsr.vcf