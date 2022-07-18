#!/bin/bash
# Get coverage stats.

# Source environment with all necessary software.
module load anaconda
source activate snakemake

# Read from command line: bam, ref genome.
ref=${1} 
bam=${2}
cov_stats=${3}

# Collect stats with Picard
picard CollectWgsMetrics \
R=${ref} \
I=${bam} \
O=${cov_stats}