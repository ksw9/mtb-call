######################
#### PGY analysis ####
######################

# Set up environmemnt
ref=/labs/jandr/walter/tb/data/refs/H37Rv.fa
DATA_DIR=/labs/jandr/walter/tb/data/
PROCESS_DIR=/labs/jandr/walter/tb/mtb/
SNPEFF_DIR=/labs/jandr/walter/tb/data/refs/snpEff/
SCRIPTS_DIR=${PROCESS_DIR}workflow/scripts/
XML_DIR=/labs/jandr/walter/tb/pgy/xml/
IQTREE_DIR=/labs/jandr/walter/tb/pgy/msa/iqtree/
module add anaconda/3_2022.05
module add IQ-TREE/2.2.0 

source activate snakemake
cd $PROCESS_DIR

# SCG items to create files for Paraguay analysis: mtb_sequencing_summary.Rmd #

#### Metadata org ####
## move metadata to cluster
rclone copy box:Box/TB/spillover/paraguay/metadata/pgy_metadata.csv /labs/jandr/walter/tb/pgy/metadata/

#### MSA ####
fasta_list=/labs/jandr/walter/tb/pgy/msa/pgy_fasta_list_single.txt
msa=/labs/jandr/walter/tb/pgy/msa/pgy_msa.fa
snps=/labs/jandr/walter/tb/pgy/msa/pgy_msa_snps.fa

# Rename fasta with duplicate name: 
dup_name=results/corrida1510/35/fasta/35_bwa_H37Rv_gatk_qfilt.fa
sed -i 's/>351/>35_corrida1510/g' $dup_name

while read line
 do #echo $line
 cat $line
done < ${fasta_list} >  ${msa}

rclone copy $msa box:Box/TB/spillover/paraguay/msa/

# Snp-sites
conda activate snp-sites
snp-sites ${msa} > ${snps} # issues with alignment length only when using the -m option. exclude this and it works fine.

# Is the problem with interleaved FASTA adding spaces? 14,100 SNP sites.
seqkit stats $snps

# IQ-Tree - find best model and fit tree
#/Applications/iqtree-1.6.12-MacOSX/bin/iqtree -s Box/Box/TB/spillover/paraguay/msa/pgy_msa.fa -st DNA -nt 8 -m HKY

# Submit to run model selection
msa=$snps
model=MFP+ASC # model search for only the models with ASC
prefix=pgy
sbatch --mem 100G -e ${prefix}.er -o ${prefix}.out --account=jandr -t 48:00:00 ${SCRIPTS_DIR}run_iqtree.sh $snps $model ${IQTREE_DIR}${prefix}

# Add boostraps - can use multiple cores (ntasks; does not matter where they are assigned)
prefix=pgy_bs
model=TVM+F+ASC+R4
sbatch --mem 100G -e ${prefix}.er -o ${prefix}.out --ntasks=16 --account=jandr -t 48:00:00 ${SCRIPTS_DIR}run_iqtree.sh $snps $model ${IQTREE_DIR}${prefix} 1000

# There is a lot of missingness in the data - try trimming alignment to see if it changes things. 
# Need to rerun with a less intensive depth filter # Also confirm that PPE filter is being applied. 

##############
#### MSA2 ####
##############
# Rename fasta with duplicate name: 
dup_name=results/corrida1510/35/fasta/35_bwa_H37Rv_gatk_qfilt.fa
sed -i 's/>35/>35_corrida1510/g' $dup_name

fasta_list=/labs/jandr/walter/tb/pgy/msa/pgy_fasta_list_single.txt
msa=/labs/jandr/walter/tb/pgy/msa/pgy_msa_080322.fa
snps=/labs/jandr/walter/tb/pgy/msa/pgy_msa_snps_080322.fa

# Update PPE masking with script. 
fasta_list=/labs/jandr/walter/tb/pgy/msa/pgy_fasta_list_single.txt
bed=/labs/jandr/walter/varcal/data/refs/ppe_gagneux_0based_snpeff.bed.gz

while read line
 do echo $line
 sbatch -A jandr -t 20 /labs/jandr/walter/tb/mtb/workflow/scripts/mask_ppe.sh $line $bed
done < ${fasta_list}

# List of masked fastas
sed 's/qfilt.fa/qfilt_masked.fa/g' $fasta_list > /labs/jandr/walter/tb/pgy/msa/pgy_fasta_list_single_masked.txt
fasta_list_masked=/labs/jandr/walter/tb/pgy/msa/pgy_fasta_list_single_masked.txt

# Create MSA.
while read line
 do #echo $line
 cat $line
done < ${fasta_list_masked} >  ${msa}

rclone copy $msa box:Box/TB/spillover/paraguay/msa/

# Snp-sites
conda activate snp-sites
snp-sites ${msa} > ${snps} # issues with alignment length only when using the -m option. exclude this and it works fine.

rclone copy $snps box:Box/TB/spillover/paraguay/msa/

# Is the problem with interleaved FASTA adding spaces? 14,100 SNP sites.
seqkit stats $snps

# IQ-Tree - find best model and fit tree
#/Applications/iqtree-1.6.12-MacOSX/bin/iqtree -s Box/Box/TB/spillover/paraguay/msa/pgy_msa.fa -st DNA -nt 8 -m HKY

cd /labs/jandr/walter/tb/pgy/msa/iqtree
# Submit to run model selection
msa=$snps
model=MFP+ASC # model search for only the models with ASC
prefix=pgy_single_masked
sbatch --mem 100G -e ${prefix}.er -o ${prefix}.out --account=jandr -t 48:00:00 ${SCRIPTS_DIR}run_iqtree.sh $snps $model ${IQTREE_DIR}${prefix}

# Wait for model selection to complete and add boostraps - can use multiple cores (ntasks; does not matter where they are assigned)
prefix=pgy_single_masked_bs
model=K3Pu+F+ASC+R5
sbatch --mem 100G -e ${prefix}.er -o ${prefix}.out --ntasks=16 --account=jandr -t 48:00:00 ${SCRIPTS_DIR}run_iqtree.sh $snps $model ${IQTREE_DIR}${prefix} 1000

# ## Instead, combine VCF files and then extract snps
# vcf_list='/labs/jandr/walter/tb/pgy/msa/pgy_vcf_list.txt'
# combined_vcf=/labs/jandr/walter/tb/pgy/msa/pgy_combined.vcf.gz
# bcftools merge --file-list ${vcf_list} -o ${combined_vcf} -O z --force-samples
# 
# # Need to reheader VCF files before merging them. 
# https://github.com/samtools/bcftools/issues/823

##################################################
#### Collate all lineage information from PGY ####
##################################################

# Copy all JSON files from tb-profiler to common directory. 
while read file; do 
  cp $file ../pgy/tb-profiler/results/
done < ../pgy/tb-profiler/tb-profiler.txt 
cd ../pgy/tb-profiler/results/
tb-profiler collate 

# Get info on percentage of ahpC mutation
cd /labs/jandr/walter/tb/pgy/tb-profiler
grep '"Rv2428","ahpC"'  ../../mtb_tgen/results/*/*/stats/*_bwa_H37Rv_lineage.csv > ahpC_profile.csv

#########################
#### Fit BEAST trees ####
#########################
conda activate snp-sites

# For cluster specific msas, select snps only. 
for msa in msa/cluster{2,4,5,11}_filt.fa ; do
  snps=${msa/.fa/_snps.fa}
  snp-sites ${msa} > ${snps}
done 

# Do the same for lineage4.4.1.1
msa=/labs/jandr/walter/tb/pgy/msa/lineage4.4.1.1_filt.fa
snps=${msa/.fa/_snps.fa}
snp-sites ${msa} > ${snps}
# Do the same for lineage4.1.2.1
msa=/labs/jandr/walter/tb/pgy/msa/lineage4.1.2.1_filt.fa
snps=${msa/.fa/_snps.fa}
snp-sites ${msa} > ${snps}
# Do the same for lineage4.3.3
msa=/labs/jandr/walter/tb/pgy/msa/lineage4.3.3_filt.fa
snps=${msa/.fa/_snps.fa}
snp-sites ${msa} > ${snps}

cd ${XML_DIR}
module load BEAST2
module load BEAGLE
module load tracer

# Move XMLs to correct directories in ${XML_DIR}filt/

# To fix clock, need to update the operators in the XML files to not sample clock rate. There are 2 operators that need to be changed.
for xml in ${XML_DIR}filt/{constant,skyline}Fixed/*xml; 
  do echo $xml; 
  # Sed to do string replacement only on lines containing a pattern.
  sed -ie '/StrictClockRateScaler/ s/weight="3.0"/weight="0.0"/g' ${xml} 
  sed -ie '/strictClockUpDownOperator/ s/weight="3.0"/weight="0.0"/g' ${xml} 
done

# Refit BEAST trees after excluding PPE genes. Running the fixed clock models again.
#for xml_path in ${XML_DIR}filt/*/*snpcor.xml ; do 
for xml_path in ${XML_DIR}filt/*/{lineage4.1.2.1,lineage4.3.3}*snpcor.xml ; do 
  cd $(dirname $xml_path)
  cd $(dirname $xml_path)
  xml=$(basename $xml_path)
  echo $xml
  sbatch -e ${xml}.out -o ${xml}.out --account=jandr -t 48:00:00 ${SCRIPTS_DIR}run_beast.sh ${xml} resume
done

# Resample trees.
for tre in ${XML_DIR}filt/*/*snps.trees; do
  echo $tre
  logcombiner --burnin 10 -renumber -resample 90000 -log $tre -o ${tre/.trees/_sampled.trees} &
done  

# For lineages, longer burn-in. 
for tree in ${XML_DIR}filt/*/{lineage4.1.2.1,lineage4.3.3}*snps.trees ; do 
  echo $tre
  logcombiner --burnin 25 -renumber -resample 90000 -log $tre -o ${tre/.trees/_sampled.trees} &
done  

# MCC trees for plotting
for tre in ${XML_DIR}filt/*/*_sampled.trees; do 
for tre in ${XML_DIR}filt/*/{lineage4.1.2.1,lineage4.3.3}*_sampled.trees; do 
  echo $tre
  treeannotator -burnin 0 -heights median ${tre} ${tre/_sampled.trees/_median_mcc.txt} &
done

# Cluster 4 skyline/fixedClock did not converge -- resumed for a total of 200 million samples.
for xml_path in ${XML_DIR}filt/skylineFixed/cluster4_filt_snps_skyline_fixed_snpcor.xml ; do 
  cd $(dirname $xml_path)
  xml=$(basename $xml_path)
  echo $xml
  sbatch -e ${xml}.out -o ${xml}.out --account=jandr -t 48:00:00 ${SCRIPTS_DIR}run_beast.sh ${xml} resume
done

# For cluster 4, redo resampling and MCC trees. 
# Resample trees.
for tre in  ${XML_DIR}filt/skylineFixed/cluster4_filt*.trees ; do
  echo $tre
  logcombiner --burnin 50 -renumber -resample 90000 -log $tre -o ${tre/.trees/_sampled.trees} &
done  

# MCC trees for plotting
for tre in ${XML_DIR}filt/skylineFixed/cluster4_filt*_sampled.trees; do 
  echo $tre
  treeannotator -burnin 0 -heights median ${tre} ${tre/_sampled.trees/_median_mcc.txt} &
done

# MultiType tree
cd /labs/jandr/walter/tb/pgy/xml/multiType

# To fix clock, need to update the operators in the XML files to not sample clock rate. There are 2 operators that need to be changed.
for xml in *xml; 
  do echo $xml; 
  # Sed to do string replacement only on lines containing a pattern.
  sed -ie '/StrictClockRateScaler/ s/weight="3.0"/weight="0.0"/g' ${xml} 
  sed -ie '/strictClockUpDownOperator/ s/weight="3.0"/weight="0.0"/g' ${xml} 
done

# Submit BEAST trees. 
for xml_path in *snpcor.xml ; do 
  cd $(dirname $xml_path)
  xml=$(basename $xml_path)
  echo $xml
  sbatch -e ${xml}.out -o ${xml}.out --account=jandr --mem 500G -t 48:00:00 ${SCRIPTS_DIR}run_beast.sh ${xml} 
done


####################
#### SRA UPLOAD ####
####################

## Copy files locally
rclone copy "/labs/jandr/walter/tb/pgy/metadata/pgy_sra_data_081222.tsv" box:Box/TB/spillover/paraguay/seq/
rclone copy "/labs/jandr/walter/tb/pgy/metadata/pgy_biosample_data_081222.tsv" box:Box/TB/spillover/paraguay/seq/

## List of files to rename
rename_list=/labs/jandr/walter/tb/pgy/metadata/file_dictionary.csv
key_dir=/tmp/key/
key_path=${key_dir}key.txt
upload_dir=/tmp/sra_upload_0822/
mkdir $upload_dir
mkdir $key_path

## Need to rename files for SRA upload so they do not include PHI. 
while read -r a b; do
  # Remove echo if satisfied by the output
  cp "$a" "$b"
done < $rename_list

## Upload
module load aspera
ascp -i ${key_path} -QT -l100m -k1 -d ${upload_dir} subasp@upload.ncbi.nlm.nih.gov:uploads/kwalter_stanford.edu_hG59bCBT

###########################
## Look at ahpC mutation ##
###########################

# Need to redo snakemake
bcftools mpileup -r  NC_000962.3:2726119-272612 trim_kr_1_bwa_H37Rv_rg.sorted.bam --fasta-ref $ref

# Bam file is consistent with TBprofiler
NC_000962.3	2726119	.	G	A,C,<*>	0	.	DP=79;I16=30,17,3,8,2177,106357,415,18435,2820,169200,660,39600,1056,25384,250,5962;QS=0.839892,0.154707,0.00540123,0;VDB=0.0657205;SGB=-0.676189;RPBZ=0.734348;MQBZ=0;MQSBZ=0;BQBZ=-1.81626;SCBZ=-0.483779;FS=0;MQ0F=0	PL	152,0,255,255,255,255,255,255,255,255

# bcftools mpileup -r  NC_000962.3:2726119  corrida1810/104/bams/tmp/104_trim_kr_1_bwa_H37Rv_rg.sorted.bam  --fasta-ref $ref | bcftools call -mv -Ob -o test/ahpC.vcf.gz
NC_000962.3     2726119 .       G       A       116.615 .       DP=79;VDB=0.0657205;SGB=-0.676189;RPBZ=0.734348;MQBZ=0;MQSBZ=0;BQBZ=-1.81626;SCBZ=-0.483779;FS=0;MQ0F=0;AC=1;AN=2;DP4=30,17,3,8;MQ=60   GT:PL   0/1:152,0,255

# VCF is not consistent though - is it overfiltered? 

# This is not working with updated tb-profiler - is it because of delly? or does it happen when i don't use fastq?
tb-profiler profile --no_delly  --bam ../mtb_tgen/results/corrida1810/104/bams/tmp/104_trim_kr_1_bwa_H37Rv_rg.sorted.bam --csv 

# Test out TB-Profiler - this works when running on fastq, 
r1=data/pgy/corrida1810/104_S4_L001_R1_001.fastq.gz 
r2=data/pgy/corrida1810/104_S4_L001_R2_001.fastq.gz 
tb-profiler profile -1 $r1 -2 $r2 -t4 --csv 

# Doesn't work on BAM because BAM is constructed with 
bam=results/corrida1810/104/bams/tmp/104_trim_kr_1_bwa_H37Rv_rg.sorted.bam
tb-profiler profile --bam $bam -t4 --csv 

# Update ref name for TB-profiler
tb-profiler update_tbdb --match_ref ${ref}

#############
#### OLD ####
#############


# Error: 8 files deleted (4 samples) in IS-1045.

IS-1045/Batch_2018-06-26/CR-654_S168_R2_001.fastq.gz
IS-1045/Batch_2018-06-26/CR-654_S22_R2_001.fastq.gz
IS-1045/Batch_2018-06-26/CR-457_S172_R2_001.fastq.gz
IS-1045/Batch_2018-06-26/CR-457_S18_R2_001.fastq.gz
IS-1045/Batch_2018-06-26/CR-512_S242_R2_001.fastq.gz
IS-1045/Batch_2018-06-26/CR-512_S33_R2_001.fastq.gz
IS-1045/Batch_2018-06-26/CR-534_S192_R2_001.fastq.gz
IS-1045/Batch_2018-06-26/CR-534_S9_R2_001.fastq.gz

# Copy files re-downloaded from Aspera (now excluding duplicates)
for f in tmp/TB-CPath-Croda-IS-1045-*; do 
  echo $f
  new_name=$(basename ${f#*TB-CPath-Croda-IS-1045-})
  cp $f IS-1045/Batch_2018-06-26/${new_name}
done

# These samples were also truncated (duplicates)
-rw-r--r--+ 1 kwalter scg_lab_jandr        52 Feb  2 15:10 T1-724-2017-30_S25_R2_001.fastq.gz
-rw-r--r--+ 1 kwalter scg_lab_jandr        52 Feb  2 15:06 T1-729-2017-54_S27_R2_001.fastq.gz
-rw-r--r--+ 1 kwalter scg_lab_jandr        53 Feb  2 14:35 T1-772-2017-814_S41_R2_001.fastq.gz
-rw-r--r--+ 1 kwalter scg_lab_jandr        52 Feb  2 14:28 T1-730-2017-47_S26_R2_001.fastq.gz
-rw-r--r--+ 1 kwalter scg_lab_jandr        53 Feb  1 19:51 T1-773-2017-812_S40_R2_001.fastq.gz
-rw-r--r--+ 1 kwalter scg_lab_jandr        53 Feb  1 19:51 T1-776-2017-868_S43_R2_001.fastq.gz
-rw-r--r--+ 1 kwalter scg_lab_jandr        53 Feb  1 19:50 T1-771-2017-825_S42_R2_001.fastq.gz

# Copy over re-downloaded from Aspera. 
for f in tmp/TB-CPath-Croda-Regina-IS-1045-*; do 
  echo $f
  new_name=$(basename ${f#*TB-CPath-Croda-Regina-IS-1045-})
  cp $f IS-1045/Batch_2018-06-26/${new_name}
done

# For duplicated samples move to duplicated data directory
mv /labs/jandr/walter/tb/data/IS-1045/Batch_2018-06-26/{CR-534_S192,CR-512_S242,CR-457_S172,CR-654_S168}_R*_001.fastq.gz tmp

########



# Also truncated by new trim_reads script

-rw-r--r--+ 1 kwalter scg_lab_jandr         51 Feb  2 15:10 IS-1018/Batch_2018-01-11/T1-724-2017-30_S9_R2_001.fastq.gz
-rw-r--r--+ 1 kwalter scg_lab_jandr         51 Feb  2 15:06 IS-1018/Batch_2018-01-11/T1-729-2017-54_S2_R2_001.fastq.gz
-rw-r--r--+ 1 kwalter scg_lab_jandr         53 Feb  2 14:35 IS-1018/Batch_2018-01-11/T1-772-2017-814_S32_R2_001.fastq.gz
-rw-r--r--+ 1 kwalter scg_lab_jandr         52 Feb  2 14:28 IS-1018/Batch_2018-01-11/T1-730-2017-47_S18_R2_001.fastq.gz
-rw-r--r--+ 1 kwalter scg_lab_jandr         53 Feb  1 19:51 IS-1018/Batch_2018-01-11/T1-773-2017-812_S41_R2_001.fastq.gz
-rw-r--r--+ 1 kwalter scg_lab_jandr         53 Feb  1 19:51 IS-1018/Batch_2018-01-11/T1-776-2017-868_S40_R2_001.fastq.gz
-rw-r--r--+ 1 kwalter scg_lab_jandr         53 Feb  1 19:50 IS-1018/Batch_2018-01-11/T1-771-2017-825_S36_R2_001.fastq.

# Copy these to IS-1018 directory. 
for f in tmp/TB-CPath-Croda-IS-1018-*; do 
  echo $f
  new_name=$(basename ${f#*TB-CPath-Croda-IS-1018-})
  echo $new_name
  cp $f IS-1018/Batch_2018-01-11/${new_name}
done

# Confirm that IS-1045 samples are not duplicated. 
