#! /bin/bash

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
export COMPUTE_SCRIPT=${SCRIPTS_DIR}/runPB_single.sh
export PDBLIST=$1
export SIZE=$2
NPROC=$3

cat $PDBLIST | xargs -n 1 -P $NPROC sh -c 'echo "### $1" && $COMPUTE_SCRIPT $1 ${SIZE}' sh
