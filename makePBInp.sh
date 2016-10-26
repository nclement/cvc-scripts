PDB=$1

POTON=$2
ENON=$3
VOLON=$4
FORON=$5

echo "computePotential ${POTON}"
echo "computeEnergy ${ENON}"
echo "computeVolumePotential ${VOLON}"
echo "computeForces ${FORON}"
echo "pqrFile PQR/${PDB}.pqr"
echo "rawnFile RAWN/${PDB}.rawn"
echo "potentialFile PBOUT/${PDB}.potential"
echo "outputPrefix PBOUT/${PDB}"
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

