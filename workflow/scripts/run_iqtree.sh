#!/bin/bash

module add IQ-TREE/2.2.0 

# Read from command line.
msa=$1
model=$2
prefix=$3
bootstraps=$4

mem_limit=100G

# Run IQ-Tree - no bootstraps if # bootstraps is not specified
if [ -z "$bootstraps" ]; then 
  iqtree2 -s ${msa} --seqtype DNA -T AUTO -m ${model} --prefix ${prefix} -mem ${mem_limit}
else
# Add bootstraps if # bootstraps is specified
  iqtree2 -s ${msa}  --seqtype DNA -T AUTO -m ${model} --prefix ${prefix} -mem ${mem_limit} -B ${bootstraps} -bnni
fi