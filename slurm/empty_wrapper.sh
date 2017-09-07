#! /bin/bash

# 1: <SUBMISSION_SCRIPT>
# 2: <JOB_NAME>
# 3: <OUTPUT_NAME>
# 4: [Job time]

SUB_SCRIPT=$1
JOB_NAME=$2
OUT_NAME=$3
JOB_TIME=${4:-24}

if [ $JOB_TIME -le 4 ]; then
  HR=4
elif [ $JOB_TIME -le 6 ]; then
  HR=6
elif [ $JOB_TIME -le 8 ]; then
  HR=8
elif [ $JOB_TIME -le 12 ]; then
  HR=12
elif [ $JOB_TIME -le 18 ]; then
  HR=18
else
  HR=24
fi

echo "Running with $HR hours"

SUBMISSION=$SUB_SCRIPT sbatch --export=SUBMISSION -J $JOB_NAME -o $OUT_NAME.o%j -N 1 -n 1 /work/01872/nclement/scripts/slurm/empty_job_${HR}h.sh
