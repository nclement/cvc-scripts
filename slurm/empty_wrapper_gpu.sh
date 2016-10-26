#! /bin/bash

# 1: <SUBMISSION_SCRIPT>
# 2: <JOB_NAME>
# 3: <OUTPUT_NAME>

sbatch -J $2 -o $3.o%j <(cat /work/01872/nclement/scripts/slurm/empty_job_gpu.sh; echo "# Everything before this line can be ignored"; cat $1)
