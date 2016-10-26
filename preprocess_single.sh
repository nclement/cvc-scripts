SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#MolSurf=/work/01872/nclement/MolSurf/release/bin/MolSurf
#PDB2PQR=/work/01872/nclement/software/pdb2pqr-src-2.1.0/pdb2pqr.py
PDB=$1

mkdir -p PQR
mkdir -p F2D
mkdir -p RAW
mkdir -p RAWN
mkdir -p TXT
mkdir -p QUAD

  
	i=$PDB
	j=PQR/$(echo "$i" | sed s/\.pdb$/.pqr/g)
	k=F2D/$(echo "$i" | sed s/\.pdb$/.f2d/g)
	l=RAW/$(echo "$i" | sed s/\.pdb$/.raw/g)
	m=TXT/$(echo "$i" | sed s/\.pdb$/.txt/g)
	n=RAWN/$(echo "$i" | sed s/\.pdb$/.rawn/g)
	o=QUAD/$(echo "$i" | sed s/\.pdb$/.quad/g)

	echo "generating pqr $i to $j"
	echo python $PDB2PQR --ff=amber "$i" "$j"
	python $PDB2PQR --ff=amber "$i" "$j"

	echo "Executing: ${MolSurf} -surfaceAtoms $j $m"
	${MolSurf} -surfaceAtoms "$j" "$m" 

	echo "Executing: ${MolSurf} -generateF2d $j $m $k 0 0"
	${MolSurf} -generateF2d "$j" "$m" "$k" 0 0

	echo "Executing: ${MolSurf} -surfaceUsingAdaptiveGrid $j $l 128"
	${MolSurf} -surfaceUsingAdaptiveGrid "$j" "$l" 128

	echo "Executing: ${MolSurf} -normals -average $l $n"
	${MolSurf} -normals -average "$l" "$n"

	echo "Executing: ${MolSurf} --aSpline -quad $n gaussian 1 $o"
	${MolSurf} -aSpline -quad "$n" gaussian 1 "$o"
