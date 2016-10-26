#! /bin/bash

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
COMPUTE_SCRIPT=${SCRIPTS_DIR}/computeAllEnergy.sh
PDB=$1
SIZE=$2
PB=$3

if [[ "$#" -eq "3" ]]; then
	$COMPUTE_SCRIPT ${PDB} ${SIZE} 1
else
	$COMPUTE_SCRIPT ${PDB} ${SIZE}
fi
