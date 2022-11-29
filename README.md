# *M. tuberculosis*  variant identification

Snakemake pipeline for *M. tuberculosis* variant identification from short-read data, mapping with bwa, variant identification with GATK. 

## Usage

Clone repository from GitHub.
```
git clone https://github.com/ksw9/mtb-call.git
```

Navigate to your project root directory. 

Create a conda environment with snakemake. Installation instructions [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
```
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Activate snakemake environment
```
source activate snakemake 
snakemake --help
```

 Update the config file so that Snakemake uses the correct sample list as input. The test sample list is config/test_data.tsv.	

 List snakemake jobs that have not yet completed, but don't run.
```
snakemake -np
```

 Run snakemake, specifying cores to use and use conda. 
```
snakemake --cores all --use-conda
```

This will stop at the Kraken taxonomic filtering step because the step for loading the Kraken database into memory doesn't work. Instead submit to the cluster. 

Submit snakemake to the cluster.
```
nohup snakemake -j 5 -k --cluster-config config/cluster_config.yaml --use-conda --rerun-triggers mtime --rerun-incomplete --cluster \
"sbatch -A {cluster.account} --mem={cluster.memory} -t {cluster.time} --cpus-per-task {threads} --error {cluster.error} --output {cluster.output} " \
> runs/snakemake_test_data.out & 
```

Monitor jobs running on the cluster.
```
sq
```

Look at entire output once jobs are completed.
```
cat runs/snakemake_test_data.out
```

Each step also produces a log, for troubleshooting a specific step. 
```
cat results/test_data/test/bams/test_bwa_H37Rv_map.log
```
 
If snakemake runs into an error or if a run is prematurely terminated, the directory will be locked.
```
snakmake --unlock
```

```

## Directory structure
Results will be organized in the results directory, in sub-directories organized by sequencing batch and then sample name.


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
│   ├── test_data/test (example organized by sequencing batch, then sample) 
|   │   ├──trim
|   │   ├──kraken
|   │   ├──bams
|   │   ├──vars
|   │   ├──fasta
|   │   ├──stats

