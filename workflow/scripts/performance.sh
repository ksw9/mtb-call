###############################
#### Pipeline performance  ####
###############################

cd /labs/jandr/walter/tb/mtb

WRK_DIR=/labs/jandr/walter/tb/mtb/working/performance/
ref=/labs/jandr/walter/tb/data/refs/H37Rv.fa

## Combine all fastas with H37Rv ref.
# Fasta files have * input at lowQual sites. Replace with N. 
for fasta in results/perform/CDC1551_clean_{1..20}/fasta/*_bwa_H37Rv_gatk*.fa; do 
  echo $fasta
  fasta_out=${fasta/.fa/_fmt.fa}
  sed 's/\*/N/g' $fa > ${fasta_out}
done
  
## Strange error with fasta files
cat ${ref} results/perform/CDC1551_clean_{1..20}/fasta/*_bwa_H37Rv_gatk_fmt.fa > ${WRK_DIR}fastas/CDC1551_bwa_H37Rv_gatk.fa
cat ${ref} results/perform/CDC1551_clean_{1..20}/fasta/*_bwa_H37Rv_gatk_PPEmask_fmt.fa > ${WRK_DIR}fastas/CDC1551_bwa_H37Rv_gatk_PPEmask.fa
cat ${ref} results/perform/CDC1551_clean_{1..20}/fasta/*_bwa_H37Rv_gatk_qPPEmask_fmt.fa > ${WRK_DIR}fastas/CDC1551_bwa_H37Rv_gatk_qPPEmask.fa

# How does bbmap work? Need to update this step--however, this extends time dramatically. 
 Allow name to be a substring of the read name (which includes /1 or /2 in these simulated paired-end reads)
NB552392:20:HWCHGBGXF:1:11101:18760:1067

# For each of the VCFs, select only variant sites to measure performance.
for replicate in {1..20}; do
  echo $replicate
  vcf=results/perform/CDC1551_clean_${replicate}/vars/CDC1551_clean_${replicate}_bwa_H37Rv_gatk.vcf.gz ; 
  output=working/performance/CDC1551_clean_${replicate}_variants.vcf.gz
  bcftools view -i 'type="SNP"' ${vcf} -O z -o ${output}
done


######################
#### Fasta tools  ####
######################

# Count nucleotides in sequence
seqtk comp

# Count nucleotides with grep 
grep -v '#' $fasta | grep -o -E 'A|C|G|T|N|-' | sort | uniq -c 

# Issue is that * are being placed into the FASTA file. 
