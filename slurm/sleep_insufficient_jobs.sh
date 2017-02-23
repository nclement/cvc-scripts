#! /bin/bash

run_jobs=`squeue -u nclement | wc -l`
# Let's keep one job open so we can submit something else!!
while [ $run_jobs -gt 49 ]; do
  run_jobs=`squeue -u nclement | wc -l`
  echo "Sleeping for a minute to prevent starting too many jobs... [$run_jobs]"
  sleep 60
done
