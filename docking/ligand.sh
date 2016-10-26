DOCKING_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR=${DOCKING_SCRIPTS_DIR}/../
source ${SCRIPTS_DIR}/Makefile.def

#MolSurf=/work/01872/nclement/MolSurf/release/bin/MolSurf
#PDB2PQR=/work/01872/nclement/software/pdb2pqr-src-2.1.0/pdb2pqr.py
FIXPQR=${SCRIPTS_DIR}/fix_pqr.pl
# The moving molecule is the ligand.

LIGAND=$1

	j=$(echo "$LIGAND" | sed s/\.pdb$/.pqr/g)
	k=$(echo "$LIGAND" | sed s/\.pdb$/.f2d/g)
	l=$(echo "$LIGAND" | sed s/\.pdb$/.raw/g)
	m=$(echo "$LIGAND" | sed s/\.pdb$/.txt/g)
	n=$(echo "$LIGAND" | sed s/\.pdb$/.rawn/g)
	o=$(echo "$LIGAND" | sed s/\.pdb$/.quad/g)

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
	echo "Executing: $MolSurf -surfaceUsingAdaptiveGrid $j $l 128"
	$MolSurf -surfaceUsingAdaptiveGrid "$j" "$l" 128
fi

if [ ! -f $n ]; then
	echo "Executing: $MolSurf -normals -average $l $n"
	$MolSurf -normals -average "$l" "$n"
fi

if [ ! -f $o ]; then
	echo "Executing: $MolSurf --aSpline -quad $n gaussian 1 $o"
	$MolSurf -aSpline -quad "$n" gaussian 1 "$o"
fi
