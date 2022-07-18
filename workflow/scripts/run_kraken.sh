#!/bin/bash
# Takes trimmed fastq, selects only reads corresponding to Mycobacterium genus or Mycobacterium tuberculosis species.
module load bbmap
# Kraken database.
DBNAME=/labs/jandr/walter/varcal/data/refs/kraken/

# Read from command line paired end 1 and paired end 2 - trimmed, zipped fastq files.
p1=$1
p2=$2
filt1=$3
filt2=$4
report=$5

# name outputs
prefix=$(basename ${p1/_trim_1.*})
LOGDIR=$(dirname $report)/

# run kraken to taxonomically classify paired-end reads and write output file.
kraken2 --db ${DBNAME} --paired --gzip-compressed --threads 8 --report ${report} --use-names ${p1} ${p2} --output ${LOGDIR}${prefix}.out

# select for reads directly assigned to Mycobacterium genus (G) (taxid 1763), reads assigned directly to Mycobacterium tuberculosis complex (G1) (taxid 77643), and reads assigned to Mycobacterium tuberculosis (S) and children. *this includes reads assigned to Mtb and those assigned to the genus, but not a different species.

grep -E 'Mycobacterium \(taxid 1763\)|Mycobacterium tuberculosis' ${LOGDIR}${prefix}.out | awk '{print $2}' > ${LOGDIR}${prefix}_reads.list

# use bbmap to select reads corresponding to taxa of interest.
filterbyname.sh  int=false in1=${p1} in2=${p2} out1=${filt1} out2=${filt2} names=${LOGDIR}${prefix}_reads.list include=true overwrite=true

# Replace bbmap with grep to filter reads according to read names.
#xzgrep -A3 -f ${LOGDIR}${prefix}_reads.list $p1 > ${filt1}
#xzgrep -A3 -f ${LOGDIR}${prefix}_reads.list $p2 > ${filt2}

# test if kraken filtered fastq files are the same length
#len1=$(cat ${filt1}  | wc -l)
#len2=$(cat ${filt2} | wc -l)  

#echo $len1 $len2

# Write error to log.
if [ $len1 != $len2 ]; then
  echo ${prefix} error: paired-end files are of different lengths. 
fi  

# Summarize Kraken statistics. 
/labs/jandr/walter/tb/scripts/kraken_stats.sh ${report} 

# Error handling
if [ "$?" != "0" ]; then
	echo "[Error]" $LINENO "Kraken ${p1} failed!" 1>&2
	exit 1
fi
