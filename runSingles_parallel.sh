#!/bin/bash

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
export COMPUTE_SCRIPT=${SCRIPTS_DIR}/computeAllEnergy.sh

# Each PDB should have the trailing ".pdb" on the end
export PDBLIST=$1
export SIZE=$2
NPROC=$3
PB=$4 #Optional, will also compute PB energy if anything is supplied here.

if [[ "$#" -eq "4" ]]; then
	# Will run PB energy as well
	cat $PDBLIST | xargs -n 1 -P $NPROC sh -c 'echo "### $1" && $COMPUTE_SCRIPT $1 ${SIZE} 1' sh
else
	cat $PDBLIST | xargs -n 1 -P $NPROC sh -c 'echo "### $1" && $COMPUTE_SCRIPT $1 ${SIZE}' sh
fi
