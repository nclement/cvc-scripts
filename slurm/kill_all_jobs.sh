#! /bin/bash

for job in $(squeue -u nclement | tr -s " " | cut -f 2 -d ' ' | grep -v JOBID); do 
  echo "Killing job $job"
  scancel $job; 
done
