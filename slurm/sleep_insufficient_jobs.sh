#! /bin/bash

run_jobs=`squeue -u nclement | wc -l`
while [ $run_jobs -gt 50 ]; do
  run_jobs=`squeue -u nclement | wc -l`
  echo "Sleeping for a minute to prevent starting too many jobs..."
  sleep 60
done
