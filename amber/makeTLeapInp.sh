PDB=$1
PDB_SHORT=`echo $PDB | sed 's/.pdb$//'`

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"/../
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def

# This is the previous stuff that doesn't actually work.
## VAC=${PDB}_vac
## echo source $AMBER_leaprc_path
## echo $PDB_SHORT = loadPdb $PDB
## echo solvateBox $PDB_SHORT WATBOX216 10 # add explicit waters
## echo saveAmberParm $PDB_SHORT AMBER/$VAC.top AMBER/$VAC.crd
## echo Quit
echo source $AMBER_leaprc_path
echo pdb = loadPdb $PDB
echo set default pbradii bondi
echo solvateBox pdb TIP3PBOX 10
echo saveAmberParm pdb AMBER/${PDB}_vac.top AMBER/${PDB}_vac.crd
echo Quit
