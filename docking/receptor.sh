DOCKING_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR=${DOCKING_SCRIPTS_DIR}/../
source ${SCRIPTS_DIR}/Makefile.def

#MolSurf=/work/01872/nclement/MolSurf/release/bin/MolSurf
#PDB2PQR=/work/01872/nclement/software/pdb2pqr-src-2.1.0/pdb2pqr.py
PQR2F2D=${SCRIPTS_DIR}/docking/f2dgen/pqrTof2d
FIXPQR=${SCRIPTS_DIR}/fix_pqr.pl

# Need to load the f2dgen module before using this
# The static molecule is the receptor.

RECEPTOR=$1

	j=$(echo "$RECEPTOR" | sed s/\.pdb$/.pqr/g)
	k=$(echo "$RECEPTOR" | sed s/\.pdb$/_1.7.f2d/g)
	l=$(echo "$RECEPTOR" | sed s/\.pdb$/.raw/g)
	m=$(echo "$RECEPTOR" | sed s/\.pdb$/_imp.raw/g)
	n=$(echo "$RECEPTOR" | sed s/\.pdb$/.rawn/g)
	o=$(echo "$RECEPTOR" | sed s/\.pdb$/.quad/g)

if [ ! -f $j ]; then
	echo "generating pqr $RECEPTOR to $j"
	python $PDB2PQR --ff=amber "$RECEPTOR" "$j"
fi

if [ ! -f $k ]; then
	echo "generating f2d $RECEPTOR to $k"
	echo $PQR2F2D "$j" 1.7
	$PQR2F2D "$j" 1.7
	#$FIXPQR $k > $k.tmp
	#cp $k.tmp $k
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
