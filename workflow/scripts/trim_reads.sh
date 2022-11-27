#!/bin/bash
# Trim paired-end fastq files for mapping. Script takes two paired-end FASTQ files, trim adapters with TrimGalore and remove the NextSeq poly-G runs with cutadapt.
# Not used

# set up environment
#source activate snakemake

# Read from command line: ref genome, fastq 1, fastq 2.
p1=$1  # fastq 1
p2=$2  # fast1 2
trim1=$3 # trimmed/adapter-removed fastq 1
trim2=$4 # trimmed/adapter-removed fastq 2
platform=${5:-miseq}

# Quit if trimmed files already exist
{
if [ -f "$trim1" ]; then
    echo "Trimmed file exists!"
    exit 0
fi
 }
 
# # Print input and output filenames.
# echo $p1 $p2 
# echo 'output files' $trim1 $trim2
# 
# # Make temp directory for intermediate files. 
# mkdir -p $(dirname $p1)/tmp
TMP_DIR=/tmp/
# 
# # Name temporary files.
# tmp1=${TMP_DIR}$(basename ${p1/.fastq.gz/_val_1.fq.gz})
# tmp2=${TMP_DIR}$(basename ${p2/.fastq.gz/_val_2.fq.gz})


#### Trim READS ####
echo "trimming reads"

# trip adapters; use default stringency (overlap of 1-bp with adapter is trimmed)
if platform == 'nextseq'; then
  trim_galore  --gzip --paired ${p1} ${p2} --nextseq 20 --output_dir ${TMP_DIR} 2>&1 
elif 
  trim_galore  --gzip --paired ${p1} ${p2} --quality 20 --output_dir ${TMP_DIR} 2>&1 
fi

#trim_galore  --gzip --stringency 3 --paired ${p1} ${p2} --output_dir ${TMP_DIR} 2>&1 
# -q Trims low-quality (Phred scaled base quality < 20) in addition to adapter removal
# --length Sets min read length at 70; bwa mem is optimized for reads >=70 bp. (Abbott)

# second, try cutadapt to remove next seq poly-G's
#cutadapt -f 'fastq' --nextseq-trim=20  --minimum-length=20 --pair-filter=any -o ${trim1} -p ${trim2} ${tmp1} ${tmp2} 2>&1 
# removes poly-G tails -NovaSeq. In those instruments, a ‘dark cycle’ (with no detected color) encodes a G. However, dark cycles also occur when when sequencing “falls off” the end of the fragment. The read then contains a run of high-quality, but incorrect `G` calls 

# Test if both paired-end reads have same number of lines.
# len1=$(zcat $trim1 | wc -l)
# len2=$(zcat $trim2 | wc -l)

# Create error if two files are different lengths.
# if [ $len1 != $len2 ]; then   
#    echo ${prefix}'_1.fq'  
#    error: paired-end files are of different lengths. 
# fi  

# remove extra files
rm ${tmp1} ${tmp2}

# Error handling
if [ "$?" != "0" ]; then
	echo "[Error]" $LINENO "Read trimming ${p1} failed!" 1>&2
	exit 1
fi