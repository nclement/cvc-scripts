#! /bin/bash
SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
# Load the modules.
source $SCRIPTS_DIR/modules_mol
#MolSurf=/work/01872/nclement/MolSurf/release/bin/MolSurf
#SCRIPTS_DIR=/work/01872/nclement/scripts
#PDB2PQR=/work/01872/nclement/software/pdb2pqr-src-2.1.0/pdb2pqr.py

PDB=$1
SIZE=$2

mkdir -p PQR
mkdir -p RAW

echo python $PDB2PQR --ff=amber $PDB PQR/${PDB}.pqr
python $PDB2PQR --ff=amber $PDB PQR/${PDB}.pqr

echo "${MolSurf} -surfaceUsingAdaptiveGrid PQR/${PDB}.pqr RAW/${PDB}.raw ${SIZE}"
${MolSurf} -surfaceUsingAdaptiveGrid PQR/${PDB}.pqr RAW/${PDB}.raw ${SIZE}

echo "${MolSurf} -removeInteriorPockets RAW/${PDB}.raw RAW/${PDB}-a.raw"
${MolSurf} -removeInteriorPockets RAW/${PDB}.raw RAW/${PDB}-a.raw

echo "${MolSurf} -qualityImprove RAW/${PDB}-a.raw RAW/${PDB}.raw"
${MolSurf} -qualityImprove RAW/${PDB}-a.raw RAW/${PDB}.raw

