PQR=$1

echo <<EOF
molPQR $PQR
epsilon 0.4
numThreads 16
numCudaDevice 1
numCudaThreads 1
methodLeaf ANLS
maxLeafAtoms 24
minRadius 2.0
threadType FLEX3
forceSixAtomType true
cudaBlockSize 128
minInterAtomDist 1.0
cilkPortion 0.3
EOF
