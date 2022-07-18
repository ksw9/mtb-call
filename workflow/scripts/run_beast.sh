#!/bin/bash
# Run BEAST

# Load BEAST and BEAGLE
module load BEAST2
module load BEAGLE

# read from command line.
xml=$1
resume_params=$2

# If no resume variable, start run afresh.
# Run BEAST
if [ -z "$resume_params" ]; then 
  echo 'starting new BEAST run'
  beast -beagle -overwrite -threads 8 $xml     
else
  # Else, resume previous run. 
  echo 'resuming run'
  beast -beagle -threads 8 -resume $xml 
fi