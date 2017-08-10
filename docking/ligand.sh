DOCKING_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR=${DOCKING_SCRIPTS_DIR}/../
HBOND_DIR=$SCRIPTS_DIR/hbond
source ${SCRIPTS_DIR}/Makefile.def

#MolSurf=/work/01872/nclement/MolSurf/release/bin/MolSurf
#PDB2PQR=/work/01872/nclement/software/pdb2pqr-src-2.1.0/pdb2pqr.py
FIXPQR=${SCRIPTS_DIR}/fix_pqr.pl
# The moving molecule is the ligand.

LIGAND=$1
QUALITY=${2:-256}

	j=$(echo "$LIGAND" | sed s/\.pdb$/.pqr/)
	k=$(echo "$LIGAND" | sed s/\.pdb$/.f2d/)
	l=$(echo "$LIGAND" | sed s/\.pdb$/.raw/)
	m=$(echo "$LIGAND" | sed s/\.pdb$/.txt/)
	n=$(echo "$LIGAND" | sed s/\.pdb$/.rawn/)
	o=$(echo "$LIGAND" | sed s/\.pdb$/.quad/)

# If we're going to be using hbond filter, need to do some other stuff.
if [ $USE_HBOND ]; then
  # Will create new files called *_pnon.{pdb,psf}
  PMIN=${LIGAND%.pdb}_pnon.pdb
  PMIN_PSF=${PMIN%.pdb}.psf
  PMIN_MOL2=${PMIN%.pdb}.mol2
  NO_TEST=true $HBOND_DIR/runHbondSingle.sh $LIGAND
  # Change the inputs to be the outputs of hbond.
  mv $PMIN $LIGAND
  mv $PMIN_PSF ${LIGAND%.pdb}.psf
  mv $PMIN_MOL2 ${LIGAND%.pdb}.mol2
else
  echo "Not using hbond filter [$USE_HBOND]"
fi

if [ ! -f $j ]; then
	echo "generating pqr $LIGAND to $j"
	python $PDB2PQR --ff=amber "$LIGAND" "$j"
fi

if [ ! -f $m ]; then
	echo "Executing: $MolSurf -surfaceAtoms $j $m"
	$MolSurf -surfaceAtoms "$j" "$m" 
fi

if [ ! -f $k ]; then
	echo "Executing: $MolSurf -generateF2d $j $m $k 0 0"
	$MolSurf -generateF2d "$j" "$m" "$k" 0 0
#	echo $FIXPQR $k \> $k.tmp
#	$FIXPQR $k > $k.tmp
#	echo mv $k.tmp $k
#	mv $k.tmp $k
fi

if [ ! -f $l ]; then
	echo "Executing: $MolSurf -surfaceUsingAdaptiveGrid $j $l $QUALITY"
	$MolSurf -surfaceUsingAdaptiveGrid "$j" "$l" $QUALITY
fi

if [ ! -f $n ]; then
	echo "Executing: $MolSurf -normals -average $l $n"
	$MolSurf -normals -average "$l" "$n"
fi

if [ ! -f $o ]; then
	echo "Executing: $MolSurf --aSpline -quad $n gaussian 1 $o"
	$MolSurf -aSpline -quad "$n" gaussian 1 "$o"
fi
