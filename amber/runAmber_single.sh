#! /bin/bash
SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
#SCRIPTS_DIR=/work/01872/nclement/scripts
source $SCRIPTS_DIR/modules_amber
#module load intel/13.0.2.146
#module load mvapich2/1.9a2
#module load amber
#module list

# Must have the correct modules loaded (amber)

PDB=$1
NCYC=${2:-1000} # If you want to specify the number of cycles, do so here.
VAC=${PDB}_vac
PDB_SHORT=`echo $PDB | sed 's/.pdb$//'`

mkdir -p AMBER
mkdir -p INP
mkdir -p ERR

echo "Processing ${PDB}..."
${SCRIPTS_DIR}/makeTLeapInp.sh ${PDB} > INP/${PDB}.tleap
# This needs to work.
if tleap -s -f INP/${PDB}.tleap | grep "FATAL" > ERR/${PDB}.tleap.error; then
  # Remove the offending atoms
  $SCRIPTS_DIR/fixAmber_errors.pl $PDB ERR/${PDB}.tleap.error $PDB
  # Try it again
  if tleap -s -f INP/${PDB}.tleap | grep "FATAL"; then
    # Otherwise, hard quit
    echo "FATAL ERROR: Could not generate tleap"
    exit
  fi
fi

${SCRIPTS_DIR}/makeAmberInp.sh $NCYC > INP/${PDB}.amberin

## echo ibrun sander -O -i INP/${PDB}.amberin -o AMBER/${PDB}.amberout -c AMBER/${VAC}.crd -p AMBER/${VAC}.top -r AMBER/${VAC}.min.crd
## # Stampede has the following function (multi-threaded)
## #ibrun pmemd.MPI -O -i INP/${PDB}.amberin -o AMBER/${PDB}.amberout -c AMBER/${VAC}.crd -p AMBER/${VAC}.top -r AMBER/${VAC}.min.crd
#echo ibrun sander.MPI -O -i INP/${PDB}.amberin -o AMBER/${PDB}.amberout -c AMBER/${VAC}.crd -p AMBER/${VAC}.top -r AMBER/${VAC}.min.crd
#ibrun sander.MPI -O -i INP/${PDB}.amberin -o AMBER/${PDB}.amberout -c AMBER/${VAC}.crd -p AMBER/${VAC}.top -r AMBER/${VAC}.min.crd

echo mpirun pmemd.MPI -O -i INP/${PDB}.amberin -o AMBER/${PDB}.amberout -c AMBER/${VAC}.crd -p AMBER/${VAC}.top -r AMBER/${VAC}.min.crd -x AMBER/${PDB}.min.nc -inf AMBER/${PDB}.min.info
mpirun pmemd.MPI -O -i INP/${PDB}.amberin -o AMBER/${PDB}.amberout -c AMBER/${VAC}.crd -p AMBER/${VAC}.top -r AMBER/${VAC}.min.crd -x AMBER/${PDB}.min.nc -inf AMBER/${PDB}.min.info

# Extract the new file.
ambpdb -p AMBER/${VAC}.top -c AMBER/${VAC}.min.crd > AMBER/${PDB_SHORT}_ambermin.pdb
# Lonestar needs the following:
# sander -O -i INP/${PDB}.amberin -o AMBER/${PDB}.amberout -c AMBER/${VAC}.crd -p AMBER/${VAC}.top -r AMBER/${VAC}.min.crd
# # Extract the new file.
# ambpdb -p AMBER/${VAC}.top -c AMBER/${VAC}.min.crd > AMBER/${PDB_SHORT}_ambermin.pdb
