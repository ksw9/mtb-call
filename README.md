# Mtb variant calling, 9/2/2022

Pipeline for *M. tuberculosis* variant identification from short-read data.

## Usage
```
# Navigate to root directory. 
# Identify files that have not been produced (no run)
snakemake -np

# Run snakemake
snakemake
```

## Directory structure

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
