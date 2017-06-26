HBOND_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR="$HBOND_SCRIPTS_DIR/../"
source $SCRIPTS_DIR/Makefile.def

PDB=$1
CHAINS=${2:-?}

	echo "Processing ${PDB}..."
  lc_PDB=$(echo $PDB | perl -ne 'print lc')
  if [ $lc_PDB != $PDB ]; then 
    echo "WARNING: PDB file should be lower case!!!"
    cp $PDB $lc_PDB
    PDB=$lc_PDB
  fi

	CPDB=$(echo "${PDB}" | sed s/\.pdb$/_cmin.pdb/g)
	PSF=$(echo "${PDB}" | sed s/\.pdb$/_cmin.psf/g)
	MOL=$(echo "${PDB}" | sed s/\.pdb$/_cmin.mol2/g)
	BASE=$(echo "${PDB}" | sed s/\.pdb$//g)
	HBOND=$(echo "${PDB}" | sed s/\.pdb$/.hbond/g)
	OUT=$(echo "${PDB}" | sed s/\.pdb$/.out/g)

	echo "perl ./pdbprep.pl ${PDB}"
	perl $HBOND_SCRIPTS_DIR/pdbprep.pl ${PDB} 

	echo "perl ./pdbchm.pl ${BASE} $CHAINS"
	perl $HBOND_SCRIPTS_DIR/pdbchm.pl ${BASE} $CHAINS --rtf=${HBOND_SCRIPTS_DIR}/pdbamino.rtf --prm=${HBOND_SCRIPTS_DIR}/parm.prm --charmm=${CHARMM}

	echo "obabel -ipdb ${CPDB} -omol2 -O ${MOL}"
	$OBABEL -ipdb ${CPDB} -omol2 -O ${MOL} > /dev/null

	#echo "testHbond ${CPDB} ${PSF} parm.prm pdbamino.rtf empty.pdb atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec ${MOL} > ${OUT}"
	echo ${F2DOCK_REFACTORED_BIN}/testHbond ${CPDB} ${HBOND} ${PSF} ${HBOND_SCRIPTS_DIR}/parm.prm ${HBOND_SCRIPTS_DIR}/pdbamino.rtf ${HBOND_SCRIPTS_DIR}/empty.pdb ${HBOND_SCRIPTS_DIR}/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec ${MOL} \> ${OUT}
	${F2DOCK_REFACTORED_BIN}/testHbond ${CPDB} ${HBOND} ${PSF} ${HBOND_SCRIPTS_DIR}/parm.prm ${HBOND_SCRIPTS_DIR}/pdbamino.rtf ${HBOND_SCRIPTS_DIR}/empty.pdb ${HBOND_SCRIPTS_DIR}/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec ${MOL}
