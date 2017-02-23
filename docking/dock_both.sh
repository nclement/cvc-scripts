
DOCKING_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
source $DOCKING_SCRIPTS_DIR/../Makefile.def

LIGAND=${DOCKING_SCRIPTS_DIR}/ligand.sh
RECEPTOR=${DOCKING_SCRIPTS_DIR}/receptor.sh
GEN_PARAM=${DOCKING_SCRIPTS_DIR}/gen-param-file.sh
#F2D=/work/01872/nclement/cvc_fresh/F2Dock-refactored/Release/bin/F2Dock
# Previous version of D2Dock gives better scoring.
F2D_prev=/work/01872/nclement/F2Dock_prev/src/runF2Dock.sh
F2D=${F2DOCK_REFACTORED}

#module load gcc/4.7.1
module load mvapich2/2.1
module load fftw3
ml

LIG=$1 # Needs to be a PDB, but if the PQR exists, just use .pdb
REC=$2 # Needs to be a PDB, and needs PDB to exist
OUT=$3 # Output file name. Should be .txt, and will add .inp to the input filename
RMSD_ATOMS=$4 # Optional, for adding RMSD. Should be ligand (moving) RMSD atoms
# Optional parameters are:
#   RES_CONT_FILTER
#      controls the residue contact filter (breaks docking with non-
#      proteins). To use this, do:
#        RES_CONT_FILTER=true dock_both.sh LIGAND RECEPTOR OUT RMSD_ATOMS
#      Values are 'true' or 'false' (or not defined)
#   TYPE
#      the type of complex to dock.
#      Default value is "U" (for Unknown). Should be E, A, or U
#   USE_PREV
#      Uses a previous version of F2Dock that will print more scores.

# Quit if the outfile already exists
has=`grep -c "END PEAKS" ${OUT}`
if [[ "$has" -eq "1" ]]; then
  exit
fi

# Set the type if it exists.
if [ -z ${TYPE+x} ]; then
  echo "No type specified in the environment TYPE=[$TYPE]"
  TYPE="U"
fi

# See if we should use F2Dock previous version.
if ! [ -z ${USE_PREV+x} ]; then
  F2D=$F2D_prev
  echo "Using previous version of F2dock"
fi

$LIGAND $LIG
$RECEPTOR $REC
echo $GEN_PARAM $LIG $REC $TYPE $OUT \> ${OUT%.txt}.inp
$GEN_PARAM $LIG $REC $TYPE $OUT > ${OUT%.txt}.inp
if [[ "$#" -eq "4" ]]; then
  echo "rmsdAtoms $RMSD_ATOMS" >> ${OUT%.txt}.inp
fi
# If this is defined, include it.
if ! [ -z ${RES_CONT_FILTER+x} ]; then
  echo applyResidueContactFilter ${RES_CONT_FILTER} >> ${OUT%.txt}.inp
fi
echo $F2D ${OUT%.txt}.inp
$F2D ${OUT%.txt}.inp
