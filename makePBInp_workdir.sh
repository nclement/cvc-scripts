WORKDIR=$1
PDB=$2
SURF=$3

POTON=$4
ENON=$5
VOLON=$6
FORON=$7

echo "computePotential ${POTON}"
echo "computeEnergy ${ENON}"
echo "computeVolumePotential ${VOLON}"
echo "computeForces ${FORON}"
echo "pqrFile ${WORKDIR}/PQR/${PDB}.pqr"
echo "rawnFile ${WORKDIR}/RAWN/${SURF}.rawn"
echo "potentialFile ${WORKDIR}/PBOUT/${PDB}.potential"
echo "outputPrefix ${WORKDIR}/PBOUT/${PDB}"
echo "ionConc 0.0"
echo "epsilonE 80"
echo "epsilonI 2"
echo "temperature 300"
echo "linearSolver GMRES"
echo "discretizationMethod Collocation"
echo "bieFormulation NBIE"
echo "geometry ASpline"
echo "quadratureLevel 1"
echo "quadOrderFarField 3"
echo "quadOrderNearField 1"
echo "quadOrderSingular 3"
echo "fmmAccuracy 4"
echo "rawivOutputResolution 128"


