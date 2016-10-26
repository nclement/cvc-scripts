
DOCKING_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
source $DOCKING_SCRIPTS_DIR/../Makefile.def

LIGAND=${DOCKING_SCRIPT_DIR}/ligand.sh
RECEPTOR=${DOCKING_SCRIPT_DIR}/receptor.sh
GEN_PARAM=${DOCKING_SCRIPT_DIR}/gen-param-file.sh
#F2D=/work/01872/nclement/cvc_fresh/F2Dock-refactored/Release/bin/F2Dock
#F2D=/work/01872/nclement/F2Dock_prev/src/runF2Dock.sh
F2D=${F2DOCK_REFACTORED}

#module load gcc/4.7.1
module load mvapich2/2.1
module load fftw3
ml

LIG=$1 # Needs to be a PDB, but if the PQR exists, just use .pdb
REC=$2 # Needs to be a PDB, and needs PDB to exist
OUT=$3 # Output file name. Should be .txt, and will add .inp to the input filename
RMSD_ATOMS=$4 # Optional, for adding RMSD. Should be ligand (moving) RMSD atoms
# An optional parameter is:
# RES_CONT_FILTER
# which will turn off the residue contact filter (breaks docking with non-
# proteins). To use this, do:
# RES_CONT_FILTER=false dock_both.sh LIGAND RECEPTOR OUT RMSD_ATOMS
# Values for this param are 'true' or 'false'

# Quit if the outfile already exists
has=`grep -c "END PEAKS" ${OUT}`
if [[ "$has" -eq "1" ]]; then
  exit
fi

$LIGAND $LIG
$RECEPTOR $REC
echo $GEN_PARAM $LIG $REC O $OUT \> ${OUT%.txt}.inp
$GEN_PARAM $LIG $REC O $OUT > ${OUT%.txt}.inp
if [[ "$#" -eq "4" ]]; then
	echo "rmsdAtoms $RMSD_ATOMS" >> ${OUT%.txt}.inp
fi
if [ -v RES_CONT_FILTER ]; then
  echo "Not applying residue contact filter!"
  echo applyResidueContactFilter ${RES_CONT_FILTER} >> ${OUT%.txt}.inp
fi
echo $F2D ${OUT%.txt}.inp
$F2D ${OUT%.txt}.inp
