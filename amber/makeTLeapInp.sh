PDB=$1
VAC=${PDB}_vac
PDB_SHORT=`echo $PDB | sed 's/.pdb$//'`

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"/../
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def

echo source ${AMBER_leaprc_path}
echo $PDB_SHORT = loadPdb $PDB
echo solvateBox $PDB_SHORT WATBOX216 10 # add explicit waters
echo saveAmberParm $PDB_SHORT AMBER/$VAC.top AMBER/$VAC.crd
echo Quit
