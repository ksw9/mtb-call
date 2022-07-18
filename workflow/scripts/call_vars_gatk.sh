#!/bin/bash
# Variant calling with GATK. Requires: (1) reference genome, (2) bam file, (3) ploidy, (4) optionally, output directory.

source activate mtb # move this to conda/envs once mamba is installed

# Read input arguments. 
ref=$1 
bam=$2
ploidy=$3 
vcf=$4

# name gVCF file
gvcf=${vcf/.vcf.gz/.g.vcf.gz}

# if no GATK dictionary, index reference genome
if [ ! -f ${ref%.*}".dict" ] ; then
echo "GATK dictionary for $ref" >&2
picard CreateSequenceDictionary \
	REFERENCE=${ref} \
	OUTPUT=${ref%.*}".dict" 
fi

# If BAM index does not exist, index BAM.
if [ ! -s ${bam}.bai ]; then 
  echo 'indexing bam' 
  samtools index ${bam}
fi
	
# Call variants with GATK 4.1, output GVCF
gatk --java-options "-Xmx100g" HaplotypeCaller \
-R ${ref} \
-ploidy ${ploidy} \
-I ${bam} \
-ERC GVCF \
-O ${gvcf}

# Index gvcf 
gatk IndexFeatureFile \
  -I ${gvcf}

# GVCF to VCF. 
gatk --java-options '-Xmx100g' GenotypeGVCFs \
-R ${ref} \
--variant ${gvcf} \
-ploidy ${ploidy} \
--include-non-variant-sites true \
--output ${vcf}
# min base quality score is 10 by default.

#Error handling
if [ "$?" != "0" ]; then
	echo "[Error]" $LINENO "GATK failed!" 1>&2
	exit 1
fi

#### PRINT OUTPUT ####
echo "Done====" >&2