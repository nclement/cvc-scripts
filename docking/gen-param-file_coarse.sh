#! /bin/bash

DOCKING_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
RPATH=${DOCKING_SCRIPTS_DIR}/Rotations
#RPATH=/work/01872/nclement/scripts/docking/Rotations
LIGAND=$1
RECEPTOR=$2
TYPE=$3 # must be "E" for enzyme, "A" for antibody, or "U" for other
OUTFILE=$4
gridSpacing=$5 # high-res is 1.2, maybe use something like 2.4?
rotFile=$6 # something like deg15.mtx for high-resolution, could use deg30.mtx for coarse
numRot=100000
numSolutions=5000 # used to be 20000
numRerank=1000 # used to be 2000

LIG_PDB=$(echo $LIGAND | sed 's/.pdb$/.f2d/')
REC_PDB=$(echo $RECEPTOR | sed 's/.pdb$/_1.7.f2d/')
PQR_L=$(echo $LIGAND | sed 's/.pdb$/.pqr/')
QUAD_L=$(echo $LIGAND | sed 's/.pdb$/.quad/')
PQR_R=$(echo $RECEPTOR | sed 's/.pdb$/.pqr/')
QUAD_R=$(echo $RECEPTOR | sed 's/.pdb$/.quad/')

echo "#"
echo "# Input files"
echo "#"
echo "staticMolecule $REC_PDB"
echo "movingMolecule $LIG_PDB"
echo "staticMoleculePQR $PQR_R"
echo "movingMoleculePQR $PQR_L"
echo "staticMoleculeQUAD $QUAD_R"
echo "movingMoleculeQUAD $QUAD_L"

echo "#"
echo "# Translational grid and FFT"
echo "#"
echo "effGridFile /work/01872/nclement/cvc_fresh/F2Dock-refactored/src/fftw-rank.txt";
echo "gridSpocing $gridSpacing"
echo "rotFile $RPATH/$rotFile"
echo "numRot $numRot"
echo "randomRotate true"

echo "#"
echo "# Output"
echo "#"
echo "numSolutions $numSolutions"
echo "outFile $OUTFILE"
    
echo "complexType $TYPE"

echo "#"
echo "# Dispersion Energy Filter"
echo "#"
echo "applyDispersionFilter false"

echo "#"
echo "# Reranking filters"
echo "#"
echo "rerank true"
echo "numRerank $numRerank"

echo "#"
echo "# Misc"
echo "#"
echo "numThreads 16"
# end gen-script.sh

