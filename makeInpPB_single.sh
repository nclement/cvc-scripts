WORKDIR=$1
PDB=$2
SURF=$3

echo "computePotential ON"
echo "computeEnergy ON"
echo "computeVolumePotential ON"
echo "computeForces OFF"
echo "pqrFile ${WORKDIR}/PQR/${PDB}.pqr"
echo "rawnFile ${WORKDIR}/RAWN/${SURF}.rawn"
echo "potentialFile ${WORKDIR}/PBOUT256/${PDB}.potential"
echo "outputPrefix ${WORKDIR}/PBOUT256/${PDB}"
echo "ionConc 0.0"
echo "epsilonE 80"
echo "epsilonI 2"
echo "temperature 300"
echo "linearSolver GMRES"
echo "discretizationMethod Collocation"
echo "bieFormulation DBIE"
echo "geometry ASpline"
echo "quadratureLevel 1"
echo "quadOrderFarField 3"
echo "quadOrderNearField 1"
echo "quadOrderSingular 3"
echo "fmmAccuracy 4"
echo "solvationVolumePotentialOnly 0"
#echo "sdfFile ${WORKDIR}/RAWIV/sdf_${PDB}.rawiv"
echo "rawivOutputResolution 128"
