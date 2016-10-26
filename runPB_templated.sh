#! /bin/bash
# Like runPB, only accepts a file in the format of runPB (only one PDB
# per line), and compares them to the input PDBCOMP. Also runs in parallel.

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
export COMPUTE_SCRIPT=${SCRIPTS_DIR}/computeAllEnergy.sh
export PDBCOMP=$1
export PDBLIST=$2
export SIZE=$3
NPROC=$4

cat $PDBLIST | xargs -n 1 -P $NPROC sh -c 'cat PQR/${PDBCOMP}.pqr PQR/${1}.pqr > PQR/${PDBCOMP}-${1}.pqr && $COMPUTE_SCRIPT ${PDBCOMP}-${1} ${SIZE}' sh
