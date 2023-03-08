# *M. tuberculosis*  variant identification--Deprecated 2023. Use: https://github.com/ksw9/mtb-vars

Snakemake pipeline for *M. tuberculosis* variant identification from short-read data, mapping with bwa, variant identification with GATK. 

## Usage

Clone repository from GitHub.
```
git clone https://github.com/ksw9/mtb-call.git
```

Navigate to your project root directory. 

Create a conda environment with snakemake. Installation instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
```
module add anaconda/3_2022.05 # On Stanford's SCG.
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Activate snakemake environment
```
source activate snakemake 
snakemake --help
```

Download SnpEff for [gene annotation](https://pcingola.github.io/SnpEff/download/).

Update the config file (config/config.yml) with the correct SnpEff path.
```
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

# Unzip file
unzip snpEff_latest_core.zip
```

Download the updated M. tuberculosis annotations. 
```
 java -jar snpEff.jar download Mycobacterium_tuberculosis_h37rv
```
 
Update the config file (config/config.yml) so that Snakemake uses the correct sample list as input. The test sample list is config/test_data.tsv.	

Create a directory to store cluster run outputs for troubleshooting.
```
mkdir logs
```

List snakemake jobs that have not yet completed, but don't run.
```
snakemake -np
```

Running snakemake locally will stop at Kraken step, due to memory requirements. 
Instead, submit to the cluster. ```--use-conda``` indicates using conda libraries specified at each step. 

```
nohup snakemake -j 5 -k --cluster-config config/cluster_config.yml --use-conda --rerun-triggers mtime --rerun-incomplete --cluster \
"sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} --error {cluster.error} --output {cluster.output} " \
> logs/snakemake_test_data.out & 
```

Monitor jobs running on the cluster.
```
sq
```

Look at entire output once jobs are completed.
```
cat run_logs/snakemake_test_data.out
```

Each step also produces a log, for troubleshooting a specific step. 
```
cat results/test_data/test/bams/test_bwa_H37Rv_map.log
```
 
If snakemake runs into an error or if a run is prematurely terminated, the directory will be locked.
```
snakmake --unlock
```

## Example data

Sampled paired-end fastq files are in the test_data directory.
An input sample .tsv file list is located at config/test_data.tsv.

## Directory structure
Results will be organized in the results directory, in sub-directories organized by sequencing batch and then sample name.

```
├── .gitignore
├── README.md
├── workflow
│   ├── rules
│   ├── envs
|   │   ├── bwa.yml
|   │   ├── gatk.yml
|   │   ├── kraken2.yml
|   │   ├── picard.yml
|   │   ├── quanttb.yml
|   │   ├── samtools.yml
|   │   ├── TBprofilerv4.3.0.yml
|   │   ├── trim-galore.yml
|   │   ├── conda-meta
│   ├── scripts
|   │   ├──process_stanford_tb.sh
|   │   ├──run_kraken.sh
|   │   ├──run_quanttb.sh
|   │   ├──reads_output.sh
|   │   ├──map_reads.sh
|   │   ├──cov_stats.sh
|   │   ├──mykrobe_predict.sh
|   │   ├──call_varsgatk.sh
|   │   ├──vcf2fasta.sh
│   ├── refs
|   │   ├──H37Rv.fa
|   │   ├──H37Rv_ppe.bed.gz
|   │   ├──ppe_hdr.txt
|   │   ├──snpEff
|   └── Snakefile
├── config
│   ├── config.yml (run/user specific parameters)
│   ├── cluster_config.yml (cluster specific parameters)
├── results
│   ├── test_data/test (example organized by sequencing batch, then sample) 
|   │   ├──trim
|   │   ├──kraken
|   │   ├──bams
|   │   ├──vars
|   │   ├──fasta
|   │   ├──stats
```
