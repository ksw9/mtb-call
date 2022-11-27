# Mtb variant calling, 9/2/2022

Pipeline for *M. tuberculosis* variant identification from short-read data.

## Usage
```
# Download scripts and test data from GitHub.

# Navigate to your project root directory. 

# Set up environment
module add anaconda

# Create a conda environment with snakemake
source activate snakemake 

# Update the config file so that Snakemake uses the correct sample list as input.

# Output snakemake jobs that need to be run
snakemake -np

# Run snakemake, specifying cores to use and use conda. 
snakemake --cores all --use-conda

# This will stop at Kraken because the step for loading the Kraken db into memory doesn't work. Instead submit to the cluster. 

# Run all samples (that are listed in samples csv).
nohup snakemake -j 5 -k --cluster-config config/cluster_config.yaml --use-conda --rerun-triggers mtime --rerun-incomplete --cluster \
"sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} --error {cluster.error} --output {cluster.output} " \
> runs/snakemake_test_data.out & 

# Monitor jobs running
sq

# Look at entire output once jobs are completed.
cat runs/snakemake_test_data.out

# Each step also produces a log, for troubleshooting a specific step
# If snakemake runs into an error or if a run is prematurely terminated, the directory will be locked
snakmake --unlock

# Identify files that have not been produced (no run)
snakemake -np

# Run snakemake
snakemake
```

## Directory structure
Results will be organized in the results directory, in sub-directories organized by sequencing batch and then sample name.

```
├── .gitignore
├── README.md
├── workflow
│   ├── rules
│   ├── envs
|   │   ├── mtb.yaml
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
|   └── Snakefile
├── config
│   ├── config.yaml
│   ├── cluster_config.yaml
├── results
│   ├── IS-1000/TB-T3.DNA.MTB-016501/ (example organized by sequencing batch, then sample) 
|   │   ├──trim
|   │   ├──kraken
|   │   ├──bams
|   │   ├──vars
|   │   ├──fasta
|   │   ├──stats
```
