#! /bin/bash
SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
#SCRIPTS_DIR=/work/01872/nclement/scripts
source $SCRIPTS_DIR/modules_amber
#module load intel/13.0.2.146
#module load mvapich2/1.9a2
#module load amber

# Must have the correct modules loaded (amber)

PDB=$1
VAC=${PDB}_vac

mkdir -p AMBER
mkdir -p INP

echo "Processing ${PDB}..."
${SCRIPTS_DIR}/makeTLeapInp.sh ${PDB} > INP/${PDB}.tleap
tleap -s -f INP/${PDB}.tleap
${SCRIPTS_DIR}/makeAmberInp.sh > INP/${PDB}.amberin
echo ibrun sander -O -i INP/${PDB}.amberin -o AMBER/${PDB}.amberout -c AMBER/${VAC}.crd -p AMBER/${VAC}.top -r AMBER/${VAC}.min.crd
ibrun sander -O -i INP/${PDB}.amberin -o AMBER/${PDB}.amberout -c AMBER/${VAC}.crd -p AMBER/${VAC}.top -r AMBER/${VAC}.min.crd
# Overwriting the previous PDB file -- is this what we want to do?
ambpdb -p AMBER/${VAC}.top < AMBER/${VAC}.min.crd > ${PDB}
