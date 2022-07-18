#!/bin/bash

# Output the number of reads and bases in each fastq file from pipeline. 

# Arguments
#samp=$1
#bat=$2
raw_fq=$1
trimmed_fq=$2
kraken_fq=$3
output=$4
#data_dir=/labs/jandr/walter/tb/data/
#results_dir=/labs/jandr/walter/tb/mtb_tgen/results/

# Specify batch
bat=$(basename $(dirname ${raw_fq}))

# Output the number of reads and total basepairs from the raw FASTQ files.
# raw_fq=$(find ${data_dir}/${bat} -name ${samp}_S*_R1_001.fastq.gz )
set -- $(pigz -dc ${raw_fq} |
    awk 'NR%4==2{c++; l+=length($0)}
          END{ print c, l }' )

raw_reads=$1
raw_bases=$2

# Output the number of reads and total basepairs from the trimmed FASTQ files.
# trimmed_fq=${results_dir}/${bat}/${samp}/trim/${samp}_trim_1.fq.gz
set -- $(pigz -dc ${trimmed_fq} |
    awk 'NR%4==2{c++; l+=length($0)}
          END{ print c, l }' )

trim_reads=$1
trim_bases=$2

# Output the number of reads and total basepairs from the Kraken filtered FASTQ files.
# kraken_fq=${results_dir}/${bat}/${samp}/kraken/${samp}_trim_kr_1.fq.gz
set -- $(pigz -dc ${kraken_fq} |
    awk 'NR%4==2{c++; l+=length($0)}
          END{ print c, l }' )

kraken_reads=$1
kraken_bases=$2

# Output to file
#samp_output=${results_dir}/${bat}/${samp}/stats/${samp}_read_counts.txt
echo 'batch,sampl,fastq_1,raw_reads,raw_bases,trim_reads,trim_bases,kraken_reads,kraken_bases' > ${output}
echo "${bat},${samp},${raw_fq},${raw_reads},${raw_bases},${trim_reads},${trim_bases},${kraken_reads},${kraken_bases}" >> ${output}

echo "${bat},${samp},${raw_fq},${raw_reads},${raw_bases},${trim_reads},${trim_bases},${kraken_reads},${kraken_bases}"