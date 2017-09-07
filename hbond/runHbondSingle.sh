#! /bin/bash

HBOND_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR="$HBOND_SCRIPTS_DIR/../"
source $SCRIPTS_DIR/Makefile.def

PDB=$1

# If you just want to generate the input files and not generate the hydrogen
# bonds, you can do so by setting NO_TEST, like:
#
#   NO_TEST=true runHbondSingle.sh myprot.pdb
#

echo "Processing ${PDB}..."

# Make sure these stay the same:
RTF=$HBOND_SCRIPTS_DIR/prms/pdbamino.rtf
PRM=${HBOND_SCRIPTS_DIR}/prms/parm.prm
ATM_PRM=${HBOND_SCRIPTS_DIR}/prms/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec

BASE=$(echo $PDB | sed s/\.pdb$//)
PNON=$(echo "${PDB}" | sed s/\.pdb$/_pnon.pdb/)
PSF=$(echo "${PDB}" | sed s/\.pdb$/_pnon.psf/)
MOL2=$(echo "${PDB}" | sed s/\.pdb$/_pnon.mol2/)
HBOND=$(echo "${PDB}" | sed s/\.pdb$/_pnon.hbond/)
OUT=$(echo "${PDB}" | sed s/\.pdb$/_pnon.out/)

# Creates ${BASE}_pnon.pdb
echo $HBOND_SCRIPTS_DIR/protein_prep/prepare.py $PDB
$HBOND_SCRIPTS_DIR/protein_prep/prepare.py $PDB

# Creates $PSF
echo $HBOND_SCRIPTS_DIR/create_psf/create_psf ${PNON} ${RTF} ${PSF}
$HBOND_SCRIPTS_DIR/create_psf/create_psf ${PNON} ${RTF} ${PSF}

# Creates .mol2
echo $OBABEL ${PNON} -O $MOL2 \&\> ${OUT}.obabel
$OBABEL ${PNON} -O $MOL2 &> ${OUT}.obabel
	
# If this hasn't been set, run the test.
if [ -z ${NO_TEST+x} ]; then
  echo ${F2DOCK_REFACTORED_BIN}/testHbond ${PNON} ${HBOND} ${PSF} ${PRM} ${RTF} ${HBOND_SCRIPTS_DIR}/empty.pdb ${ATM_PRM} ${MOL2} \> ${OUT}
  ${F2DOCK_REFACTORED_BIN}/testHbond ${PNON} ${HBOND} ${PSF} ${PRM} ${RTF} ${HBOND_SCRIPTS_DIR}/empty.pdb ${ATM_PRM} ${MOL2} > ${OUT}
fi
