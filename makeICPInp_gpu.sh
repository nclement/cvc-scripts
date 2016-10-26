PQR=$1

echo <<EOL
molPQR $PQR
epsilon 0.5 
printStatus true
numThreads 16
numCudaDevice 1
numCudaThreads 1
methodLeaf ANLS
maxLeafAtoms 24
minRadius 1.4
threadType CUDA
forceSixAtomType true
cudaBlockSize 128
minInterAtomDist 1.0
EOL
