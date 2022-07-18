########################################
##### Process Vitoria Mtb seq data #####
########################################


#### Generate preliminary data for within-host Mtb study ####

DATA_DIR=/labs/jandr/walter/tb/data/vitoria/
SCRIPTS_DIR=/labs/jandr/walter/tb/mtb/workflow/scripts/
PROCESS_DIR=/labs/jandr/walter/tb/vitoria/
module load anaconda; source activate snakemake

## 1. Download Colangeli et al. 2020 data. SNA:PRJNA475130. (Pair identifier is Library name.)
cd $DATA_DIR
while read sra; do 
  echo $sra
  sbatch ${SCRIPTS_DIR}download_sra.sh ${sra} ${DATA_DIR}
done < PRJNA475130_accessions.txt

## 2. Run snp calling pipeline. 
cd ${PROCESS_DIR}

# Run Snakemake on the cluster. 
today=$(date +"%Y-%m-%d")
nohup snakemake -j 500 -k --cluster-config config/config.yaml --cluster "sbatch -A {cluster.account} --mem={cluster.mem} -t {cluster.time} --cpus-per-task {threads} -o {cluster.output}" > runs/snakemake_${today}.out & 

## 3. Test if household pairs share more within host iSNVs than epidemiologically unlinked pairs. 