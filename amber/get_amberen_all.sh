#! /bin/bash
#SBATCH -n 1                 # total number of mpi tasks requested
#SBATCH -p normal            # queue (partition) -- normal, development, etc.
#SBATCH -t 24:00:00          # run time (hh:mm:ss)
#SBATCH --mail-user=nathanlclement@gmail.com
#SBATCH --mail-type=begin    # email me when the job starts
#SBATCH --mail-type=end      # email me when the job finishes
#SBATCH -A P-MRI2

# This will get the amber-style energy of all files in the regex.
#
# Run it with something like:
#  get_amberen_all.sh '/work/01872/nclement/energy_test/zdock/1CGI_r_?.pdb'

pdbregex=$1

for fullfile in $pdbregex; do 
  file=`basename $fullfile`
  # Copy it here
  cp $fullfile $file
  echo "XX Running $file"; 
  # Get rid of het atoms
  sed -i -e "/^HETATM/d" $file
  # This will make sure it doesn't fail.
  /work/01872/nclement/scripts/amber/runAmber_single.sh $file;
  # Get the energy from this file.
  echo $file $( $WORK/scripts/amber/getAmberEnergy.sh $file )
  # Then delete the original file.
  rm $file
done
