#! /bin/bash

# Except ignore any idev jobs.
for job in $(squeue -u nclement | grep -v "development" | tr -s " " | cut -f 2 -d ' ' | grep -v JOBID); do 
  echo "Killing job $job"
  scancel $job; 
done
