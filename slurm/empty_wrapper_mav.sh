#! /bin/bash

# 1: <SUBMISSION_SCRIPT>
# 2: <JOB_NAME>
# 3: <OUTPUT_NAME>

SUBMISSION=$1 sbatch --export=SUBMISSION -J $2 -o $3.o%j <(cat /work/01872/nclement/scripts/slurm/empty_job_mav.sh; echo "# Everything before this line can be ignored"; cat $1)
