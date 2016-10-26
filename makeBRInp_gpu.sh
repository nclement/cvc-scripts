PQR=$1
QUAD=$2

echo <<EOF
molPQR $PQR
molQUAD $QUAD
epsilon 0.5 
printStatus true
numThreads 16
numCudaDevice 1
numCudaThreads 1
methodLeaf ANLS
maxLeafAtoms 28
minRadius 1.3
threadType CILK
maxBornRadius 1000.0
forceSixAtomType true
writeBornRadius true
cudaBlockSize 128
maxCudaKernelLoop 32
maxCudaThreadBlock 2048
EOF
