#########################################
##### Process Stanford Mtb seq data #####
#########################################

#### Set up work-space
qlogin -m 100 -c 12

# Set up environmemnt
ref=/labs/jandr/walter/tb/data/refs/H37Rv.fa
DATA_DIR=/labs/jandr/walter/tb/data/
PROCESS_DIR=/labs/jandr/walter/tb/mtb/
SNPEFF_DIR=/labs/jandr/walter/tb/data/refs/snpEff/
SCRIPTS_DIR=/labs/jandr/walter/tb/mtb/workflow/scripts/
module add anaconda/3
module add IQ-TREE/2.2.0 
conda activate snakemake

## Push to github
git add -A 
git commit -m 'adding scripts'
git push/git pull

######################################
#### Download data from BaseSpace ####
######################################

cd ${DATA_DIR}
module load aws-cli

# List completed runs
$HOME/bin/bs list project

# Download runs
$HOME/bin/bs download project --name=MT01_MtB_Baits-2021-10-01rerun2 -o MT01_MtB_Baits-2021-10-01rerun2
$HOME/bin/bs download project --name=MT01_MtB_Baits-2021-10-06rerun2 -o MT01_MtB_Baits-2021-10-06rerun2
$HOME/bin/bs download project --name=MT01_MtB_Baits-2021-09-17 -o MT01_MtB_Baits-2021-09-17   
$HOME/bin/bs download project --name=MT02_MTB_2021-10-29 -o MT02_MTB_2021-10-29
$HOME/bin/bs download project --name=MT03_MTB_2021-11-24 -o MT03_MTB_2021-11-24 
$HOME/bin/bs download project --name=MT04_A10-12_MTB_2021-12-10 -o MT04_A10-12_MTB_2021-12-10
$HOME/bin/bs download project --name=MT06-2022-03-04_withEnvCult -o MT06-2022-03-04_withEnvCult

# List samples within a single project
$HOME/bin/bs list biosample --project-name=MT01_MtB_Baits-2021-10-06rerun2 

# Move FASTQs out of subdirectories for easier processing
for dir in *; do 
 echo $dir
 mv ${dir}/*/*fastq.gz ${dir}/
done

############################
#### Download  PGY data ####
############################
cd /labs/jandr/walter/tb/data/pgy

rclone copy dropbox:corrida1510 corrida1510
rclone copy dropbox:'CORRIDA 1509' corrida1509
rclone copy dropbox:'CORRIDA 1810' corrida1810
rclone copy dropbox:'CORRIDA 1811' corrida1811
rclone copy dropbox:'CORRIDA 2210' corrida2210
rclone copy dropbox:'CORRIDA 2309' corrida2309
rclone copy dropbox:'CORRIDA 3010' corrida3010
rclone copy dropbox:'CORRIDA 0710' corrida0710
rclone copy dropbox:'CEDIC 3011' cedic3011

# Copy shared google drive files 
rclone copy drive:Paraguay pgy1

# Copy additional sequences from Natalie Weiler
rclone copy drive:Secuencias weiler 
rclone copy drive:"Segunda carpeta de Secuencias" weiler2

# corrida1509 and pgy1 are duplicates. move pgy to duplicated directory
mv pgy1 duplicated/

#####################################
#### Download data from SeqMatik ####
#####################################
SEQMATIK_DIR=SU022
# List data available for run1
aws s3 ls --profile wasabi --endpoint-url=https://s3.wasabisys.com s3://seqmatic-data-releases/katharine_walter/20220206_SU022_Andrews_MT04Next_IndAB/

# Download to data directory
aws s3 sync --profile wasabi --endpoint-url=https://s3.wasabisys.com s3://seqmatic-data-releases/katharine_walter/20220206_SU022_Andrews_MT04Next_IndAB/ ${SEQMATIK_DIR}

SEQMATIK_DIR=SU025
# List data available for run1
aws s3 ls --profile wasabi --endpoint-url=https://s3.wasabisys.com s3://seqmatic-data-releases/katharine_walter/20220421_SU025_MT05A/

# Download to data directory
aws s3 sync --profile wasabi --endpoint-url=https://s3.wasabisys.com s3://seqmatic-data-releases/katharine_walter/20220421_SU025_MT05A/ ${SEQMATIK_DIR}

SEQMATIK_DIR=SU027
# List data available for run1
aws s3 ls --profile wasabi --endpoint-url=https://s3.wasabisys.com s3://seqmatic-data-releases/katharine_walter/20220629_SU027_MT05B/

# Download to data directory
aws s3 sync --profile wasabi --endpoint-url=https://s3.wasabisys.com s3://seqmatic-data-releases/katharine_walter/20220629_SU027_MT05B/ ${SEQMATIK_DIR}

##############################
#### Process capture runs #### old (run in mtb directory). updated Snakemake file outputs: mtb_tgen
##############################
cd ${PROCESS_DIR}

# Run Snakemake on the cluster. 
today=$(date +"%Y-%m-%d")
nohup snakemake -j 500 -k --cluster-config config/cluster_config.yaml  \
--cluster "sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} " > runs/snakemake_${today}.out &

# Look at run errors. 
runs/snakemake_${today}.out

nohup snakemake $file -j 500 -k --cluster-config config/config.yaml --cluster \
 'sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} --cpus-per-task {threads} --output "runs/slurm/slurm-%j.out"' > runs/snakemake_${today}.out &
 
# Try dry-run locally: 
snakemake -np process/vars/10561-Room_S8_bwa_MTB_ancestor_reference_gatk_vqsrfilt.vcf.gz

# Try on cluster for single sample
nohup snakemake results/8328-BC_S9/vars/8328-BC_S9_bwa_H37Rv_gatk_ann.vcf.gz -j 1 -k --cluster-config config/cluster_config.yaml  --cluster "sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} " > runs/snakemake_${today}.out &

# Test run. 
snakemake results/MT01_MtB_Baits-2021-10-06rerun2/4514-lighter-384/stats/4514-lighter-384_read_counts.txt  -j 1 -k --cluster-config config/cluster_config.yaml --cluster "sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} " > runs/snakemake_${today}.out &


# snakemake question
# https://stackoverflow.com/questions/69578275/snakemake-combine-inputs-over-one-wildcard

# Dry run of a single rule. 
snakemake -R per_samp_stats -np

# Combine coverage thus far. 
sed -n '1'p  $(ls process/stats/*combined_stats.csv | head -n1) > combined_cov_stats.csv
for file in *process/stats/*combined_stats.csv ; 
  do sed -n '2'p $file 
done >> combined_cov_stats.csv

# Copy coverage stats locally. 
rclone copy combined_cov_stats.csv box:Box/TB/seq_data/cov/ 

# Buid DAG of jobs *need to comment out print command; this causes issue with dot command.
snakemake --dag  results/fasta/4514-*fa | dot -Tsvg > process/plots/snakemake_dag.svg

###############################
#### Environmental samples ####
###############################
cd ${PROCESS_DIR}

## Create MSA of 7 environmental samples + 1 culture sample
cat $ref results/{4514,10561,8328}*/fasta/*qfilt.fa > enviro/enviro_msa.fa

## Get snps fasta and VCF
snp-sites -m enviro/enviro_msa.fa -o enviro/enviro_msa_snps.fa
snp-sites -v enviro/enviro_msa.fa -o enviro/enviro.vcf

# Replace chromosome name and zip
cat enviro/enviro.vcf | awk '{gsub(/^1/,"NC_000962.3"); print}' | bgzip > enviro/enviro.vcf.gz
tabix -f enviro/enviro.vcf.gz

# Get sites of fixed differences
bcftools query -f 'NC_000962.3\t%POS\n' enviro/enviro.vcf.gz > enviro/fixed_diffs.txt

# For all environmental VCFs & H37Rv, get het sites and sites included in the SNP positions, output as iSNV VCF. 
for vcf in results/{4514,10561,8328,H37RV_S1}*/vars/*qfilt.vcf.gz ; do 
  vcf_filt=enviro/per_samp/$(basename ${vcf/.vcf.gz})_isnv.vcf.gz
  vcf_fixed_diffs=enviro/per_samp/$(basename ${vcf/.vcf.gz})_fixed.vcf.gz
  #bcftools filter --include 'AD[0:0] >= 5 & AD[0:1] >= 5 & TYPE != "indel"' $vcf | bgzip > ${vcf_filt}
  # Also get a vcf with fixed differences for each VCF. 
  bcftools view -T 'enviro/fixed_diffs.txt' ${vcf} | bgzip > ${vcf_fixed_diffs}
done

# Read this into R. /labs/jandr/walter/tb/mtb/scripts/plot_coverage.Rmd

# Collate all tb-profiler runs (directory needs to be named results)
cp results/*/stats/results/*json results/stats/results/
cd results/stats/
tb-profiler collate 

# Also collect json files from IS-1000 run
cp  ../mtb_tgen/results/IS-1000/*/stats/results/*json results/stats/results/
cd results/stats/
tb-profiler collate 

# Copy environmental plots locally. 
for file in process/plots/* ; do 
  rclone copy $file  box:Box/TB/Environmental/plots/
done

## Filter enviro VCFs to exclude MAF >5% and < 95%
for vcf in results/{4514,10561,8328,H37RV_S1}*/vars/*qfilt.vcf.gz ; do 
 echo $vcf
 vcf_filt=enviro/per_samp/$(basename ${vcf/.vcf.gz})_nohets.vcf.gz
# bcftools filter --exclude 'AD[0:0]/FORMAT/DP > 0.1 & AD[0:0]/FORMAT/DP < 0.9' --set-GTs . $vcf | bgzip > ${vcf_filt}
 bcftools filter -e 'AD[0:0]/(FORMAT/DP) > .05 && AD[0:1]/(FORMAT/DP) > 0.05' --set-GTs . $vcf | bgzip > ${vcf_filt}

done

# VCF to fasta and get pairwise differences. 
bed='/labs/jandr/walter/varcal/data/refs/ppe_gagneux_0based_snpeff.bed.gz'
for vcf in enviro/per_samp/*_nohets.vcf.gz; do
  tabix $vcf
  sample_name=$(basename ${vcf/_bwa_H37Rv_gatk_qfilt_fixed.vcf.gz})
  echo $sample_name
  fasta=enviro/per_samp/${sample_name}_fixed.fa
  ${SCRIPTS_DIR}vcf2fasta.sh $ref $vcf $sample_name $bed $fasta
done

## Get snps fasta and VCF
for samp in enviro/per_samp/*_bwa_H37Rv_gatk_qfilt_nohets.vcf.gz_fixed.fa; do
  mv $samp ${samp/_bwa_H37Rv_gatk_qfilt_nohets.vcf.gz_fixed.fa/_nohets.fa}
done

cat $ref enviro/per_samp/*nohets.fa > enviro/enviro_fixed_msa.fa
snp-sites -m enviro/enviro_fixed_msa.fa -o enviro/enviro_fixed_msa_snps.fa

## Compare to unmasked fastas.
cat results/{4514,10561,8328,H37RV_S1}*/fasta/*qfilt.fa $ref > enviro/enviro_msa.fa
snp-sites -m enviro/enviro_msa.fa -o enviro/enviro_msa_snps.fa

# Look at variants in H37Rv. 
cat $ref ../results/H37RV_S1/fasta/H37RV_S1_bwa_H37Rv_gatk_qfilt.fa > pos_control.fa
snp-sites pos_control.fa 

## Create MSA of 7 environmental samples + 1 culture sample. Add culture sample of additional Rif+ case.
cat $ref results/{4514,10561,8328}*/fasta/*qfilt.fa ../mtb_tgen/results/Stanford/11080-7651-envculture/fasta/11080-7651-envculture_bwa_H37Rv_gatk_qfilt.fa \
> enviro/enviro_culture_msa.fa

## Get snps fasta and VCF
snp-sites -m enviro/enviro_culture_msa.fa -o enviro/enviro_culture_msa_snps.fa

######################
#### Enviro Trees ####
######################
#### Enviro Trees ####
cd /labs/jandr/walter/tb/mtb/

# A. Create an MSA with all environmental sequences. ## Use updated RAXML-ng!
cat enviro/enviro_fixed_msa.fa ../mtb_tgen/results/IS-1000/*/fasta/*qfilt.fa > enviro/stanford_msa.fa 
snp-sites -m enviro/stanford_msa.fa -o enviro/stanford_msa_snps.fa

#### Updating RaxML jobs here ####
cd /labs/jandr/walter/tb/mtb/enviro/
ref=/labs/jandr/walter/varcal/data/refs/H37Rv.fa
msa=stanford_msa_snps.fa
prefix=trees/stanford_is-1000_workers
threads=8
pin=thread-nopin
brlen=mod
model=GTR
/labs/jandr/walter/tb/scripts/raxml_ng_parallel.sh \
$msa $ref $model $prefix $threads $pin $brlen
# Submit job. 
sbatch --account=jandr -t 48:00:00 /labs/jandr/walter/tb/scripts/raxml_ng_parallel.sh \
$msa $ref $model $prefix $threads $pin $brlen -o trees/stanford_is-1000.out \
 -e trees/stanford_is-1000.er
 

# Testing now with increased # of workers.
# the "local" run is the fastest. 

# Extract reads that map to rpoB point mutation. 
bam=results/10561-Food_S7/bams/10561-Food_S7_bwa_H37Rv.merged.rmdup.bam

samtools view -h ${bam} NC_000962.3:761155 | samtools fasta - # Errors with IO

# files are blocked by large tarball:
/labs/jandr/walter/tb/samps.tar.tgz

########################
#### Enviro upload ####
#######################
cd $PROCESS_DIR
module load aspera
# Copy to SRA
key=/labs/jandr/walter/lambda/metadata/aspera.openssh
file_path=/tmp/kwalter/sra_upload
mkdir -p $file_path

# Copy BAM files of interest to filepath
cp results/{10561,4514,8328}*/bams/*rmdup.bam ${file_path}

# Rename incorrectly named sample.
cd ${file_path}
mv 8328-Food_S10_bwa_H37Rv.merged.rmdup.bam 4514-Food_S10_bwa_H37Rv.merged.rmdup.bam
ascp -i ${key} -QT -l100m -k1 -d ${file_path} subasp@upload.ncbi.nlm.nih.gov:uploads/kwalter_stanford.edu_hG59bCBT

###########################
#### Process TGen runs ####
###########################
cd /labs/jandr/walter/tb/mtb_tgen
DATA_DIR=/labs/jandr/walter/tb/data/
today=$(date +'%Y_%m_%d')

sample_list=config/tgen_samples.tsv
echo 'sample,fastq_1,fastq_2','batch','run'> ${sample_list}

for file in ${DATA_DIR}IS-1000/Walter-TB/*R1_001.fastq.gz; do
  fq1=${file}
  fq2=${file/R1_001.fastq.gz/R2_001.fastq.gz}
  samp=$(basename $file)
  samp=${samp/_S*_R1_001.fastq.gz}
  dir=$(dirname $file)
  run=${dir##*/}
  bat=$(dirname $(echo ${dir/"${DATA_DIR}"}))
  echo ${samp},${fq1},${fq2},${bat},${run} >> ${sample_list}
done

# Batch 1018 - need to process -S files separately to make sure there is no duplication.
sample_list=config/tgen_samples2.tsv
echo 'sample,fastq_1,fastq_2','batch','run'> ${sample_list}

for file in ${DATA_DIR}IS-1018/*/*R1_001.fastq.gz; do
  fq1=${file}
  fq2=${file/R1_001.fastq.gz/R2_001.fastq.gz}
  samp=$(basename $file)
  samp=${samp/_S*_R1_001.fastq.gz}
  dir=$(dirname $file)
  run=${dir##*/}
  bat=$(dirname $(echo ${dir/"${DATA_DIR}"}))
  echo ${samp},${fq1},${fq2},${bat},${run} >> ${sample_list}
done

# Other batches from TGEN
sample_list=config/tgen_samples3.tsv
echo 'sample,fastq_1,fastq_2','batch','run'> ${sample_list}

for file in ${DATA_DIR}IS-{1045,1062}/*/*R1_001.fastq.gz; do
  fq1=${file}
  fq2=${file/R1_001.fastq.gz/R2_001.fastq.gz}
  samp=$(basename $file)
  samp=${samp/_S*R1_001.fastq.gz} # need to keep suffix here because of repeated samples in the same batch. 
  dir=$(dirname $file)
  run=${dir##*/}
  bat=$(dirname $(echo ${dir/"${DATA_DIR}"}))
  echo ${samp},${fq1},${fq2},${bat},${run} >> ${sample_list}
done

##################
#### SeqMatik ####
##################
cd /labs/jandr/walter/tb/mtb_tgen
DATA_DIR=/labs/jandr/walter/tb/data/

# Batch1 Seqmatk 
sample_list=config/SU022.tsv
echo 'sample,fastq_1,fastq_2','batch','run'> ${sample_list}

for file in /labs/jandr/walter/tb/data/SU022/*R1_001.fastq.gz; do
  fq1=${file}
  fq2=${file/R1_001.fastq.gz/R2_001.fastq.gz}
  samp=$(basename $file)
  samp=${samp/_S*_R1_001.fastq.gz}
  dir=$(dirname $file)
  run=${dir##*/}
  #bat=$(dirname $(echo ${dir/"${DATA_DIR}"}))
  bat=$run
  echo ${samp},${fq1},${fq2},${bat},${run} >> ${sample_list}
done


# Batch3 Seqmatk 
sample_list=config/SU025.tsv
echo 'sample,fastq_1,fastq_2','batch','run'> ${sample_list}

for file in /labs/jandr/walter/tb/data/SU025/*R1_001.fastq.gz; do
  fq1=${file}
  fq2=${file/R1_001.fastq.gz/R2_001.fastq.gz}
  samp=$(basename $file)
  samp=${samp/_S*_R1_001.fastq.gz}
  dir=$(dirname $file)
  run=${dir##*/}
  #bat=$(dirname $(echo ${dir/"${DATA_DIR}"}))
  bat=$run
  echo ${samp},${fq1},${fq2},${bat},${run} >> ${sample_list}
done

#############
#### PGY ####
#############
cd /labs/jandr/walter/tb/mtb_tgen
DATA_DIR=/labs/jandr/walter/tb/data/

# Which sample doesn't have paired end read? 

# PGY
sample_list=config/pgy.tsv # 147 samples total including both weiler directories.
echo 'sample,fastq_1,fastq_2','batch','run'> ${sample_list}

for file in ${DATA_DIR}pgy/{c,w}*/*R1_001.fastq.gz ; do
  fq1=${file}
  fq2=${file/R1_001.fastq.gz/R2_001.fastq.gz}
  
  # Only include if fq2 exists.
	if [ -f "$fq2" ]
	then
	  samp=$(basename $file)
	  samp=${samp/_S*_R1_001.fastq.gz}
   	  dir=$(dirname $file)
   	  runname=${dir##*/} # Use runname so that no run-sample is repeated. 
   	  bat=$(dirname $(echo ${dir/"${DATA_DIR}"}))
   	#bat=$run
  	  echo ${samp},${fq1},${fq2},${runname},${runname} >> ${sample_list}
	else
	  echo "${fq2} doesn't exist"
	fi
done

#  514_S14_L001_R2_001.fastq.gz  is missing an R1 file. (pgy/weiler)
# corrida2309/D_S16_L001_R1_001.fastq.gz is truncated. 
# Processing 143 smaples. 

# Get list of other pgy files with batch. 

# Move 'corrida1509' samples tmp - these are duplicates of the Google Drive file (pgy1)
mv pgy/corrida1509/* pgy/corrida1509/tmp/

# Move Google drive files to corrida1509
mv pgy/pgy1/* pgy/corrida1509/

###############################
#### Process Stanford runs ####
###############################
cd /labs/jandr/walter/tb/mtb_tgen
DATA_DIR=/labs/jandr/walter/tb/data/

# List samples
sample_list=config/MT06-2022-03-04_withEnvCult.tsv
echo 'sample,fastq_1,fastq_2','batch','run'> ${sample_list}

for file in /labs/jandr/walter/tb/data/MT06-2022-03-04_withEnvCult/*R1_001.fastq.gz ; do
  fq1=${file}
  fq2=${file/R1_001.fastq.gz/R2_001.fastq.gz}
  
  # Only include if fq2 exists.
	if [ -f "$fq2" ]; then
		samp=$(basename $file)
		samp=${samp/_S*_R1_001.fastq.gz}
   		dir=$(dirname $file)
   		bat=$(basename $dir)
		run=${dir##*/}
  		echo ${samp},${fq1},${fq2},${bat},${run} >> ${sample_list}
	else
		echo "${fq2} doesn't exist"
	fi
done

## Process earlier Stanford runs. ##

for batch in $(ls -d ${DATA_DIR}/MT0{1,2,3,4}*) ; 
 do echo ${batch}
 bat=$(basename ${batch})

	sample_list=config/${bat}.tsv
	echo 'sample,fastq_1,fastq_2','batch','run'> ${sample_list}

	for file in /labs/jandr/walter/tb/data/${bat}/*R1_001.fastq.gz ; do
	  fq1=${file}
	  fq2=${file/R1_001.fastq.gz/R2_001.fastq.gz}
  
	  # Only include if fq2 exists.
		if [ -f "$fq2" ]; then
			samp=$(basename $file)
			samp=${samp/_S*_R1_001.fastq.gz}
			dir=$(dirname $file)
			run=${dir##*/}
			#bat=$(dirname $(echo ${dir/"${DATA_DIR}"}))
			#bat=$run
			echo ${samp},${fq1},${fq2},${bat},${run} >> ${sample_list}
		else
			echo "${fq2} doesn't exist"
		fi
	done
done

# Combine sample lists of previous Stanford runs. 
sample_list='config/Stanford_MT01-04.tsv'
echo 'sample,fastq_1,fastq_2','batch','run'> ${sample_list}
cat config/MT0{1,2,3,4}*.tsv | grep -v 'sample' >> ${sample_list}

# Run Snakemake on the cluster. 
# Need to activate snakemake to run snakemake; 
source activate snakemake

today=$(date +"%Y-%m-%d") 

# Test a single sample: --use-conda
nohup snakemake -j 1 -k --cluster-config config/cluster_config.yaml --use-conda --cluster \
"sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} --error {cluster.error} --output {cluster.output} " \
results/IS-1018/T1-XX-2015-439/fasta/T1-XX-2015-439_bwa_H37Rv_gatk_qfilt.fa> runs/snakemake_${today}.out & 

# Test for updated lineage assignment & also quanttb. 



nohup snakemake -j 1 -k --cluster-config config/cluster_config.yaml --use-conda --cluster \
"sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} --error {cluster.error} --output {cluster.output} " \
results/IS-1018/T1-174-2011-715/vars/T1-174-2011-715_bwa_H37Rv_gatk_rename.vcf.gz > runs/snakemake_${today}.out & 

# Run all samples (that are listed in samples csv)
nohup snakemake -j 500 -k --cluster-config config/cluster_config.yaml --rerun-triggers mtime --cluster \
"sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} --error {cluster.error} --output {cluster.output} " \
> runs/snakemake_${today}.out & 

# Don't use the use-conda flag now -- let Jody Phelan know/ also need mamba to create the environments on the fly. Also needs mamba to be installed. 
nohup snakemake -j 500 -k --cluster-config config/cluster_config.yaml  --rerun-triggers mtime --cluster \
"sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} --error {cluster.error} --output {cluster.output} " > runs/snakemake_${today}.out & 

#########################################
#### Create 2022 sequence dictionary ####
#########################################

# Output list of files in same format as previous sequence dictionary.827 new sequences.
dir=data
seq_date=2022
country=both
new_seq_list=metadata/seq_dictionary_2022.csv
echo 'path,dir','batch','seq_date','samp','strain','seq_name','record_id','country'> ${new_seq_list}

for batch_name in data/{MT02,MT03,MT04,MT06,SU022,SU025}*/; do
  batch=$(basename $batch_name)
  echo ${batch}
  for path in ${batch_name}*R1_001.fastq.gz; do 
    echo $path
    samp=$(basename ${path/_L001_R1_001.fastq.gz})
	seq_name=${samp/_S*}
	strain=${samp/#*_}
  echo ${path},${dir},${batch},${seq_date},${samp},${strain},${seq_name},${record_id},${country}
  done
done >> ${new_seq_list}

rclone copy ${new_seq_list} box:Box/TB/spillover/metadata/dna/ 


# Output list of files in same format as previous sequence dictionary. 127 additional samples. 
dir=data
seq_date=2022
country=Paraguay
new_seq_list=pgy/metadata/pgy_2022.csv
echo 'path,dir','batch','seq_date','samp','strain','seq_name','record_id','country'> ${new_seq_list}

for batch_name in data/pgy/{c,w}*; do
  batch=$(basename $batch_name)
  #echo ${batch}
  for path in ${batch_name}/*R1_001.fastq.gz; do 
    #echo $path
    samp=$(basename ${path/_L001_R1_001.fastq.gz})
	seq_name=${samp/_S*}
	strain=${samp/#*_}
  echo ${path},${dir},${batch},${seq_date},${samp},${strain},${seq_name},${record_id},${country}
  done
done >> ${new_seq_list}

rclone copy ${new_seq_list} box:Box/TB/spillover/paraguay/metadata/
rclone copy ${new_seq_list} box:Box/TB/spillover/metadata/dna/

#########################################
#### Collate all lineage information ####
#########################################

# Copy all JSON files from tb-profiler to common directory. 

cp results/*/stats/results/*json results/stats/results/
cd results/stats/
tb-profiler collate 

#########################################
#### Collate all stats information ####
#########################################
cat results/*/*/stats/*_bwa_H37Rv_all_stats.csv > working/combined_all_stats_0522.csv

#########################################
#### Where does loss of reads occur? ####
#########################################
while read line; do
  samp=$(echo ${line} | cut -f1)
  #fastq_1=$(echo $line | cut -f2)
  #bat=$(echo $line | cut -f4) 
  echo "$samp $bat"
done < 'config/Stanford_MT01-04.tsv'

# Loop over samples and get reads output.
for bat in $(ls -d /labs/jandr/walter/tb/mtb_tgen/results/* ); do 
  for samp in $(ls ${bat}); do 
    echo $samp
    bat_name=$(basename $bat)
    sbatch -A jandr -t 5 ${SCRIPTS_DIR}reads_output.sh ${samp} ${bat_name}
  done
done

# Update Snakemake file to include this step

# Add step to combine all Snakemake outputs into a single file. 

##############
#### COV #####
##############

## Combine all coverage thus far. 
cat results/*/*/stats/*combined_stats.csv > working/combined_stats.csv

########################
#### QUANTTB Setup #####
########################
# Create a conda environment for running QUANTTB. Requires Python 2.7 
conda env export -n quanttb > envs/quanttb.yaml
conda install -c bioconda seqtk -n mtb

######################
#### snpEff set-up ####
#######################

module load snpeff

snpEff="/scg/apps/software/snpeff/4.3t/snpEff/snpEff.jar"
bed='/labs/jandr/walter/varcal/data/refs/ppe_gagneux_0based_fmt.bed'
snpeff_bed='/labs/jandr/walter/varcal/data/refs/ppe_gagneux_0based_snpeff.bed'
gff=/labs/jandr/walter/tb/data/refs/H37Rv.gff.gz

# Download TB database locally : java -jar -Xmx50g ${snpEff} download Mycobacterium_tuberculosis_h37rv -dataDir /labs/jandr/walter/tb/data/refs/snpEff/
# config file is updated with snpEff datadir.
java -jar -Xmx50g ${snpEff} download Mycobacterium_tuberculosis_h37rv -dataDir /labs/jandr/walter/tb/data/refs/snpEff/

# Create a header line for adding PPE annotation
echo -e '##FORMAT=<ID=PPE,Number=1,Type=String,Description="Located within PE/PPE genes">' > /labs/jandr/walter/tb/data/refs/ppe_hdr.txt

# Is the bed file annoting correctly? Confirmed bed file is 0/1-based.  
# Chromosome      33581   33794   Rv0031 # From Bed file

### From gff (1-based)
# zcat $gff | grep Rv0031
# Chromosome	ena	gene	33582	33794	.	+	.	ID=gene:Rv0031;biotype=protein_coding;description=Possible remnant of a transposase;gene_id=Rv0031;logic_name=ena
# Chromosome	ena	mRNA	33582	33794	.	+	.	ID=transcript:CCP42753;Parent=gene:Rv0031;biotype=protein_coding;transcript_id=CCP42753
#     
