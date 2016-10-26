SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
# Load the modules.
source $SCRIPTS_DIR/modules_mol
#source /work/01872/nclement/scripts/.module
#MolSurf=/work/01872/nclement/MolSurf/release/bin/MolSurf
#PDB2PQR=/work/01872/nclement/software/pdb2pqr-src-2.1.0/pdb2pqr.py

PDB=$1
SIZE=${2:-128}
BASE=`echo $PDB | sed 's/\.pdb$//'`

echo python $PDB2PQR --ff=amber $BASE.pdb $BASE.pqr
python $PDB2PQR --ff=amber $BASE.pdb $BASE.pqr

echo "${MolSurf} -surfaceUsingAdaptiveGrid ${BASE}.pqr ${BASE}.rawn ${SIZE}"
${MolSurf} -surfaceUsingAdaptiveGrid ${BASE}.pqr ${BASE}.rawn ${SIZE}

echo "${MolSurf} -removeInteriorPockets ${BASE}.rawn ${BASE}-a.rawn"
${MolSurf} -removeInteriorPockets ${BASE}.rawn ${BASE}-a.rawn

echo "${MolSurf} -qualityImprove ${BASE}-a.rawn ${BASE}.rawn"
${MolSurf} -qualityImprove ${BASE}-a.rawn ${BASE}.rawn

echo "${MolSurf} -aSpline -quad ${BASE}.rawn gaussian 1 ${BASE}.quad"
${MolSurf} -aSpline -quad ${BASE}.rawn gaussian 1 ${BASE}.quad
