#!/bin/bash
# Map fastq files to reference genome. Requires (1) ref genome, (2) mapping algorithm, (3) full path to read file 1, (4) full path to read file 2 (if exists), (5) prefix for output files (if paired end reads).

# Set up environment
module load picard
#source activate snakemake

# Read from command line: ref genome, fastq 1, fastq 2.
ref=${1} # 1st input is full path to reference genome
mapper=$2 # 2nd input is mapping algorithm 
p1=$3  # 3th input is full path to read 1
p2=$4  # 4th input is full path to read 2 (if paired-end)
bam=$5 # bam file.

# Define output directory. 
BAMS_DIR=$(dirname ${bam})

# Define temp directory.
TMP_DIR=${BAMS_DIR}/tmp/
mkdir -p ${TMP_DIR}

# Define prefix.
prefix=$(basename ${p1%%.*})

# Set names of intermediate files.
ref_index=${ref%.*}
ref_prefix=$(basename $ref_index)

# remove suffix
sam=${TMP_DIR}${prefix}_${mapper}_${ref_prefix}'.sam' 
rawbam=${TMP_DIR}${prefix}_${mapper}_${ref_prefix}'.bam' 
rgbam=${TMP_DIR}${prefix}_${mapper}_${ref_prefix}'_rg.bam' 
sortbam=${TMP_DIR}${prefix}_${mapper}_${ref_prefix}'_rg.sorted.bam'
echo $sam

# help message if no inputs are entered (-z checks if 1st arugment string is null).
if [ -z ${1} ]
then
  echo "This script will use Bowtie 2 to align two paired-end fastq files"
  echo "$(basename $0) <refGenome> <mapper> <dataset> <readFile1> <readFile2>"
  echo "First input is the Mtb ref genome"
  echo "Second input is the read mapper"
  echo "Third input is the dataset"
  echo "Forth input is location of file containing read 1"
  echo "Fifth input is location of file containing read 2 (if paired-end)" 
  exit
fi

# if no p2 given 
if [ -z ${p2} ]
then
  echo "read pair not specified"
fi

# get machine id_lane from SRA read identifier (old Illumina fastq format)
# if gzipped fastq
if [[ ${p1} == *.gz ]]; then
  seqid=$(zcat ${p1} | head -n 1)
else
# if not gzipped
  seqid=$(cat ${p1} | head -n 1)
fi
seqid="$(echo $seqid | cut -d' ' -f1)"
seqid="$(echo $seqid | cut -d':' -f3)"

id_lane=${seqid:-readgroup1} # default is "readgroup1"

#### MAPPING ####
# bowtie2 mapping
if [ $mapper == 'bowtie' ] || [ $mapper == 'bowtie2' ] ; then

  # if no indexing, index reference genome
  if [ ! -f ${ref%.*}".1.bt2" ] ; then
  echo "bowtie2 indexing $ref" >&2
  ${BOWTIE2_BUILD} ${ref} ${ref_index}
  fi 

  # map
  #if paired-end reads
  echo "mapping with bowtie2" >&2
  if [ ! -z ${p2} ]; then 
  bowtie2 --threads 7 -X 1100 -x ${ref_index} -1 ${p1} -2 ${p2} -S ${sam}
  # -x basename of reference index 
  # --un-gz gzips sam output
  # p is number of threads
  # end-to-end mapping is the default mode
  # -X 2000 increase maximum fragment length from default (500) to allow longer fragments (paired-end only) -X 2000 **used for roetzer data ***
  # ART simulations X = 1100 (mean fragment length = 650bp + 3 x 150-bp stdev)
  
  # if single-end reads
  elif [ -z ${p2} ]; then
    echo "single reads"
  bowtie2 --threads 7 -X 1100 -x ${ref_index} -U ${p1} -S ${sam}
  # -U for unpaired reads
  fi
  # Error handling
  if [ "$?" != "0" ]; then
    echo "[Error]" $LINENO "bowtie mapping failed ${p1}!" 1>&2
    exit 1
  fi
fi

# bwa mapping
if [ $mapper == 'bwa' ]; then
  # if no indexing, index reference genome
  if [ ! -f ${ref}".sa" ] ; then
  echo "bwa indexing $ref" >&2
  bwa index ${ref} 
  fi

  # map
  echo "mapping with bwa" >&2
  # if paired-end reads
  if [ ! -z ${p2} ]; then 
    bwa mem -t 7 ${ref} ${p1} ${p2} > ${sam}
   # -t no. of threads.
  # if single-end reads
  elif [ -z ${p2} ]; then
      echo "single reads"
    bwa mem -t 7 ${ref} ${p1} >  ${sam}
  fi
  # Error handling
  if [ "$?" != "0" ]; then
    echo "[Error]" $LINENO "bwa mapping failed ${p1}!" 1>&2
    exit 1
  fi
fi

### POST-PROCESSING ####

# Convert sam to bam 
sambamba view -t 7 -S -h ${sam} -f bam -o ${rawbam}
# -S auto-detects input format, -h includes header, -o directs output

# Add/replace read groups for post-processing with GATK
picard AddOrReplaceReadGroups \
  INPUT=${rawbam} \
  OUTPUT=${rgbam} \
  RGID=${id_lane} \
  RGLB=library1 \
  RGPU=${id_lane} \
  RGPL="illumina" \
  RGSM=${prefix}

# Sort the BAM 
sambamba sort ${rgbam}

# Index BAM
sambamba index ${sortbam}

# Remove duplicates. (-r = remove)
sambamba markdup -r -p -t7 ${sortbam} ${bam}

# Remove intermediate files.
rm ${sam} ${rawbam} ${rgbam} ${sortbam} ${sortbam}.bai

# Error handling
if [ "$?" != "0" ]; then
    echo "[Error]" $LINENO "Remove dups failed ${p1}!" 1>&2
    exit 1
fi

### PRINT OUTPUT ####
echo "Done====" >&2
