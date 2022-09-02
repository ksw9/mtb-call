#!/bin/bash
# Get coverage stats.

# Read from command line: bam, ref genome.
ref=${1} 
bam=${2}
cov_stats=${3}

# Collect stats with Picard
picard CollectWgsMetrics \
R=${ref} \
I=${bam} \
O=${cov_stats}