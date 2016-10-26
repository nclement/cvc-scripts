SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def

#MolEnergy=/work/01872/nclement/MolEnergy/release/bin/MolEnergy
#MolSurf=/work/01872/nclement/MolSurf/release/bin/MolSurf
#PDB2PQR=/work/01872/nclement/software/pdb2pqr-src-2.1.0/pdb2pqr.py
#SCRIPTS_DIR=/work/01872/nclement/scripts

PDB=$1
SIZE=$2 #256
WORKDIR=`pwd`

mkdir -p INP
mkdir -p PBOUT256
mkdir -p RAW
mkdir -p RAWN

  PDB=`echo ${PDB} | sed "s/.pdb$//"`
  if [ ! -f RAWN/${PDB}.rawn ]; then
  	echo "Doing preprocessing steps for ${PDB}..."
		echo python $PDB2PQR --ff=amber "${PDB}.pdb" PQR/"${PDB}.pqr"
		python $PDB2PQR --ff=amber "${PDB}.pdb" "PQR/${PDB}.pqr"

		echo "Executing: ${MolSurf} -surfaceUsingAdaptiveGrid PQR/${PDB}.pqr RAW/${PDB}.raw 128"
		${MolSurf} -surfaceUsingAdaptiveGrid "PQR/${PDB}.pqr" "RAW/${PDB}.raw" 128

		echo "Executing: ${MolSurf} -normals -average RAW/${PDB}.raw RAWN/${PDB}.rawn"
		${MolSurf} -normals -average "RAW/${PDB}.raw" "RAWN/${PDB}.rawn"
  fi

  echo "Processing ${PDB}..."
  ${SCRIPTS_DIR}/makeInpPB_single.sh ${WORKDIR} ${PDB} ${PDB} >${WORKDIR}/INP/${PDB}.inp
  echo ${MolEnergy} -PB ${WORKDIR}/INP/${PDB}.inp
  ${MolEnergy} -PB ${WORKDIR}/INP/${PDB}.inp
