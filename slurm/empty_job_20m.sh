#! /bin/bash -x

#SBATCH -n 1                 # total number of mpi tasks requested
#SBATCH -p normal            # queue (partition) -- normal, development, etc.
#SBATCH -t 0:20:00          # run time (hh:mm:ss)
#SBATCH --mail-user=nathanlclement@gmail.com
#SBATCH --mail-type=begin    # email me when the job starts
#SBATCH --mail-type=end      # email me when the job finishes
#SBATCH -A P-MRI2

bash $SUBMISSION
