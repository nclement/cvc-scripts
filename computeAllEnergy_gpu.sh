SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
# Load the modules.
source $SCRIPTS_DIR/modules_mol
#MolSurf=/work/01872/nclement/MolSurf/release/bin/MolSurf
#MolEnergy=/work/01872/nclement/MolEnergy/release/bin/MolEnergy
#GPUEnergy=/work/01872/nclement/OPARef
#SCRIPTS_DIR=/work/01872/nclement/scripts
#PDB2PQR=/work/01872/nclement/software/pdb2pqr-src-2.1.0/pdb2pqr.py

# Note: PDB must have the trailing ".pdb" on the end
PDB=$1
SIZE=$2
PB=$3

mkdir -p OUT
mkdir -p INP
mkdir -p PQR
mkdir -p RAWN
mkdir -p RAW
mkdir -p QUAD

if [ ! -f PQR/${PDB}.pqr ]; then
	echo python $PDB2PQR --ff=amber $PDB PQR/${PDB}.pqr
	python $PDB2PQR --ff=amber $PDB PQR/${PDB}.pqr
fi

echo "${MolSurf} -surfaceUsingAdaptiveGrid PQR/${PDB}.pqr RAWN/${PDB}.rawn ${SIZE}"
${MolSurf} -surfaceUsingAdaptiveGrid PQR/${PDB}.pqr RAWN/${PDB}.rawn ${SIZE}

echo "${MolSurf} -removeInteriorPockets RAWN/${PDB}.rawn RAWN/${PDB}-a.rawn"
${MolSurf} -removeInteriorPockets RAWN/${PDB}.rawn RAWN/${PDB}-a.rawn

echo "${MolSurf} -qualityImprove RAWN/${PDB}-a.rawn RAWN/${PDB}.rawn"
${MolSurf} -qualityImprove RAWN/${PDB}-a.rawn RAWN/${PDB}.rawn

echo "Creating RAW file"
cut -f 1,2,3 -d ' ' RAWN/${PDB}.rawn > RAW/${PDB}.raw

echo "${MolSurf} -area RAWN/${PDB}.rawn >OUT/${PDB}.area"
${MolSurf} -area RAWN/${PDB}.rawn >OUT/${PDB}.area

echo "${MolSurf} -volume RAWN/${PDB}.rawn OUT/${PDB}.volume"
${MolSurf} -volume RAWN/${PDB}.rawn OUT/${PDB}.volume

echo "${MolSurf} -aSpline -quad RAWN/${PDB}.rawn gaussian 1 QUAD/${PDB}.quad"
${MolSurf} -aSpline -quad RAWN/${PDB}.rawn gaussian 1 QUAD/${PDB}.quad

echo "${SCRIPTS_DIR}/SH/makeILJInp.sh ${PDB} >INP/${PDB}-ilj.inp"
${SCRIPTS_DIR}/makeILJInp.sh ${PDB} >INP/${PDB}-ilj.inp

echo "${MolEnergy} -internalLJ INP/${PDB}-ilj.inp >OUT/${PDB}-ilj.out"
${MolEnergy} -internalLJ INP/${PDB}-ilj.inp >OUT/${PDB}-ilj.out

echo "${MolEnergy} -internalCP INP/${PDB}-ilj.inp >OUT/${PDB}-icp.out"
${MolEnergy} -internalCP INP/${PDB}-ilj.inp >OUT/${PDB}-icp.out

echo "${SCRIPTS_DIR}/SH/makeGpolInp.sh ${PDB} >INP/${PDB}-gpol.inp"
${SCRIPTS_DIR}/makeGpolInp.sh ${PDB} >INP/${PDB}-gpol.inp

echo "${MolEnergy} -gpol INP/${PDB}-gpol.inp >OUT/${PDB}-gpol.out"
${MolEnergy} -gpol INP/${PDB}-gpol.inp >OUT/${PDB}-gpol.out

echo "${MolEnergy} -dispersion INP/${PDB}-gpol.inp >OUT/${PDB}-disp.out"
${MolEnergy} -dispersion INP/${PDB}-gpol.inp >OUT/${PDB}-disp.out



##./SH/makePBInp.sh ${PDB} >INP/${PDB}-PB.inp ON OFF OFF OFF
##${MolEnergy} -PB INP/${PDB}-PB.inp >OUT/${PDB}-pb-pot.out

##./SH/makePBInp.sh ${PDB} >INP/${PDB}-PB.inp OFF ON OFF OFF
##${MolEnergy} -PB INP/${PDB}-PB.inp >OUT/${PDB}-pb-en.out

##./SH/makePBInp.sh ${PDB} >INP/${PDB}-PB.inp OFF OFF ON OFF
##${MolEnergy} -PB INP/${PDB}-PB.inp >OUT/${PDB}-pb-vol.out

# use this if, besides computing energy, you also want the volume potentials and surfaces colored by potentials (used for visualization purposes)
if [[ "$#" -eq "3" ]]; then
	mkdir -p PBOUT256
	WORKDIR=`pwd`
	echo "Processing PB for ${PDB}..."
	${SCRIPTS_DIR}/makeInpPB_single.sh ${WORKDIR} ${PDB} ${PDB} >${WORKDIR}/INP/${PDB}.inp
	echo ${MolEnergy} -PB ${WORKDIR}/INP/${PDB}.inp
	${MolEnergy} -PB ${WORKDIR}/INP/${PDB}.inp

	# Hack for currect work--just create the energy
	#mkdir -p PBOUT
	#${SCRIPTS_DIR}/makePBInp_workdir.sh ${WORKDIR} ${PDB} ${PDB} OFF ON OFF OFF >INP/${PDB}-PB.inp
	#echo ${MolEnergy} -PB INP/${PDB}-PB.inp \>OUT/${PDB}-pb-en.out
	#${MolEnergy} -PB INP/${PDB}-PB.inp >OUT/${PDB}-pb-en.out
fi

# Clean up temporary files, just keep the RAW files
#rm RAWN/${PDB}.rawn
rm RAWN/${PDB}-a.rawn
#rm QUAD/${PDB}.quad
#rm INP/${PDB}-*.inp
