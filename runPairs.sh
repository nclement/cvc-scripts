#! /bin/bash
SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
#SCRIPTS_DIR=/work/01872/nclement/scripts/

export COMPUTE_SCRIPT=${SCRIPTS_DIR}/computeAllEnergy.sh
export PDBLIST=$1
export SIZE=$2
export NPROC=$3
PB=$4

if [[ "$#" -eq "4" ]]; then
	cat $PDBLIST | xargs -n 2 -P $NPROC sh -c 'cat PQR/${1}.pqr PQR/${2}.pqr > PQR/${1}-${2}.pqr && $COMPUTE_SCRIPT ${1}-${2} ${SIZE} 1' sh
else
	cat $PDBLIST | xargs -n 2 -P $NPROC sh -c 'cat PQR/${1}.pqr PQR/${2}.pqr > PQR/${1}-${2}.pqr && $COMPUTE_SCRIPT ${1}-${2} ${SIZE}' sh
fi
