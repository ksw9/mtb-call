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

#module add anaconda/3_2022.05
#module add IQ-TREE/2.2.0 
conda activate snakemake

## Push to github for full pipeline.
git add ... specify scripts
git commit -m 'adding scripts'
git push 

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

# Batch4 Seqmatk 
sample_list=config/SU027.tsv
echo 'sample,fastq_1,fastq_2','batch','run'> ${sample_list}

for file in /labs/jandr/walter/tb/data/SU027/*R1_001.fastq.gz; do
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

# Move 'corrida1509' samples tmp - these are duplicates of the Google Drive file (pgy1)
mv pgy/corrida1509/* pgy/corrida1509/tmp/

# Move Google drive files to corrida1509
mv pgy/pgy1/* pgy/corrida1509/


#####################
#### Validation ####
####################
cd /labs/jandr/walter/tb/mtb
DATA_DIR=/labs/jandr/walter/varcal/single/

# Which sample doesn't have paired end read? 

# Performance files
sample_list=config/perform.tsv # 147 samples total including both weiler directories.
echo 'sample,fastq_1,fastq_2','batch','run'> ${sample_list}

for file in ${DATA_DIR}CDC1551_clean_{1..20}/data/*_1.fq.gz; do
  fq1=${file}
  fq2=${file/_1.fq.gz/_2.fq.gz}
  
  # Only include if fq2 exists.
	if [ -f "$fq2" ]
	then
	  samp=$(basename $file)
	  samp=${samp/_1.fq.gz}
   	  dir=$(dirname $file)
   	  runname=perform # Use runname so that no run-sample is repeated. 
   	  bat=$(dirname $(echo ${dir/"${DATA_DIR}"}))
   	#bat=$run
  	  echo ${samp},${fq1},${fq2},${runname},${runname} >> ${sample_list}
	else
	  echo "${fq2} doesn't exist"
	fi
done


###############################
#### Process Stanford runs ####
###############################
cd /labs/jandr/walter/tb/mtb
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

######################
#### PROCESS RUNS ####
######################
# Run Snakemake on the cluster. 
# Need to activate snakemake to run snakemake; 
cd /labs/jandr/walter/tb/mtb
module add anaconda/3_2022.05
source activate snakemake

today=$(date +"%Y-%m-%d") 

# See what runs have not completed
snakemake -np --rerun-triggers mtime 

# Run all samples (that are listed in samples csv)
nohup snakemake -j 500 -k --cluster-config config/cluster_config.yaml --use-conda --rerun-triggers mtime --rerun-incomplete --cluster \
"sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} --error {cluster.error} --output {cluster.output} " \
> runs/snakemake_${today}.out & 

# Test a single sample: --use-conda (so that rule-specific conda environments are activated)
nohup snakemake -j 1 -k --cluster-config config/cluster_config.yaml --use-conda --cluster \
"sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} --error {cluster.error} --output {cluster.output} " \
results/perform/CDC1551_clean_1/kraken/CDC1551_clean_1_kr_1.fq.gz > runs/snakemake_${today}.out & 

# STUCK WITH BBMAP STEP in run_kraken.sh #### HERE ####
filterbyname.sh in=$p1 out=test.fq names=AE000516.2-2935798/1  include=true overwrite=true substring=f, ths=t
filterbyname.sh in=$p1 out=test.fq names=${LOGDIR}${prefix}_reads.list include=true overwrite=true substring=t tws=t

# Error here: 
#snakemake  results/IS-1018/T1-105-2010-139/stats/T1-105-2010-139_bwa_H37Rv_combined_stats.csv --use-conda --cores all

# CreateCondaEnvironmentException:
# Could not create conda environment from /oak/stanford/scg/lab_jandr/walter/tb/mtb_tgen/workflow/envs/mtb.yaml:
# Command:
# mamba env create --quiet --file "/oak/stanford/scg/lab_jandr/walter/tb/mtb_tgen/.snakemake/conda/5b3e765eb8210c23d169553fd7853814.yaml" --prefix "/oak/stanford/scg/lab_jandr/walter/tb/mtb_tgen/.snakemake/conda/5b3e765eb8210c23d169553fd7853814"
# Output:
# Preparing transaction: ...working... done
# Verifying transaction: ...working... failed
# 
# PaddingError: Placeholder of length '80' too short in package /oak/stanford/scg/lab_jandr/walter/tb/mtb_tgen/.snakemake/conda/5b3e765eb8210c23d169553fd7853814/bin/Rscript.
# The package must be rebuilt with conda-build > 2.0.

# Issue has previously been found here: https://github.com/ENCODE-DCC/atac-seq-pipeline/issues/103
workflow/envs/mtb.yaml

#1. Try creating a simple picard only environment. Same error. 

#2. Try changing location according to this: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#specifying-a-location-for-an-environment

# Currently, the snakemake file points to existing conda environments, not the .yaml file. 
conda env create -f envs/mtb.yaml # currently need to create the environments before running the pipeline. 

# Improper masking

## known issue with bcftools masking: https://github.com/samtools/bcftools/issues/1592

# Script to temporarily mask. Then go back to update and test Snakemake script.
bedtools maskfasta -fi $fa -bed $bed -fo masked.fa

#########################################
#### Create 2022 sequence dictionary ####
#########################################
seq_list=working/seq_dictionary_1002022.csv

# Create sequence dictionary with seq metadata about all files
cat config/Stanford_MT01-04.tsv \
  <(sed '1d' config/MT06-2022-03-04_withEnvCult.tsv ) \
  <(sed '1d' config/pgy.tsv) \
  <(sed '1d' config/SU022.tsv) \
  <(sed '1d' config/SU025.tsv) \
  <(sed '1d' config/SU027.tsv) \
  <(sed '1d' config/tgen_samples.tsv) \
  <(sed '1d' config/tgen_samples2.tsv) \
  <(sed '1d' config/tgen_samples3_nodups.tsv) > ${seq_list} # 2739 isolates

# Add allstats.csv files. Keep header only for the first file.  
sed -e '2,${/^batch/d' -e '}' results/*/*/stats/*_bwa_H37Rv_allstats.csv >  working/allstats_100422.csv 
sed -e '2,${/^batch/d' -e '}' results/*/*/stats/*_bwa_H37Rv_combstats.csv >  working/combstats_101022.csv  # Without quanttb (because quanttb has multiple lines for mixed samples)

#########################################
#### Collate all lineage information ####
#########################################
cd /labs/jandr/walter/tb/mtb
source activate TBprofiler2

# Copy all JSON files from tb-profiler to results directories; however some samples have duplicate names. Keep batch structure.  
for bat in results/IS-1062/ ; 
  do echo $bat 
    batch_name=$(basename $bat)
    mkdir working/tbprofile/${batch_name}
    cp results/${batch_name}/*/stats/results/*json working/tbprofile/${batch_name}
done

# For each batch, collate tb-profiler information.
for bat in working/tbprofile/IS-1062 ; do
  batch_name=$(basename $bat)
  echo $batch_name
  tb-profiler collate --dir working/tbprofile/${batch_name} --prefix working/tbprofile_collate/${batch_name} \
--all_variants --mark_missing 
done

# Combine all TB-profiler. Use only *txt output (not additional output files--moved to outputs dir)
awk 'BEGIN{FS="\t"; OFS=FS} {print FILENAME,$0}' working/tbprofile_collate/* > working/tbprofiler_combined_100622.csv

# For each batch, also combine quanttb results. 
awk 'BEGIN{OFS=","} {print FILENAME,$0}' results/*/*/stats/*quanttb.csv > working/quanttb_combined_100622.csv


# # File cleaning: 
# # If sample directory has been incorrectly created, move to temp directory. (IS-1018 and IS-1045)
# for samp in results/IS-1062/*; do 
#   count=`ls -1 ${samp}/stats/*allstats.csv 2>/dev/null | wc -l `
#   if [ $count == 0 ]
#   then 
#   echo "File does not exist" $samp
#   mv $samp tmp3/
#   fi 
# done | wc -l




# Output list of files in same format as previous sequence dictionary.827 new sequences.
dir=data
seq_date=2022
country=both
new_seq_list=metadata/seq_dictionary_1002022.csv
echo 'path,dir','batch','seq_date','samp','strain','seq_name','record_id','country'> ${new_seq_list}

for batch_name in data/{MT02,MT03,MT04,MT06,SU022,SU025,SU027}*/; do
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
#### Collate all stats information ####
#########################################
tail -n1 -q results/*/*/stats/*_bwa_H37Rv_all_stats.csv > working/combined_all_stats_090622.csv

cat <(head -n1 results/IS-1018/T1-XX-2015-978/stats/T1-XX-2015-978_bwa_H37Rv_all_stats.csv) working/combined_all_stats_090622.csv > working/combined_all_stats_hdr_090622.csv
 
 
cat results/*/*/stats/*_bwa_H37Rv_cov_stats.txt > working/combined_cov_stats_090222.csv
cat results/*/*/stats/*_read_counts.csv > working/combined_cov_stats_090222.csv

# Loop over samples and get reads output.
for bat in $(ls -d /labs/jandr/walter/tb/mtb/results/* ); do 
  for samp in $(ls ${bat}); do 
    echo $samp
    bat_name=$(basename $bat)
    sbatch -A jandr -t 5 ${SCRIPTS_DIR}reads_output.sh ${samp} ${bat_name}
  done
done
    paste <(sed -n '7,8'p $cov_stats )
    
    
     > {output.combined_stats}
    
    # Combine with reads output
    paste  -d ',' <( cat {output.combined_stats} | tr "\\t" "," ) {input.reads_stats} > {output.all_stats}

cov # doesn't include the sample name!



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

# Update Snakemake file to include this step -done.

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

#### Testing ####
 
## Look for character occurrence in file
grep -v '>' results/IS-1000/TB-T3.DNA.MTB-016583/fasta/TB-T3.DNA.MTB-016583_bwa_H37Rv_gatk.fa | sed -e 's/\n//g' | grep -m1 -b -o "*"


#############################
#### Performance testing ####
#############################


########################
#### Install MtbSeq ####
########################
module add anaconda/3_2022.05
# https://github.com/ngs-fzb/MTBseq_source
mamba create -n mtbseq -c bioconda mtbseq
conda activate mtbseq

gatk3-register repos/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2 


