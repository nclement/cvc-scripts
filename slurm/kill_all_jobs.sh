#! /bin/bash

for job in $(squeue -u nclement | tr -s " " | cut -f 2 -d ' '); do 
  scancel $job; 
done