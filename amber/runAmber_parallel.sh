#! /bin/bash
export AMBER_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"

export PDBLIST=$1
export NPROC=$2
export NCYC=${3:-1000}  # Number of cycles, default to 1000

cat $PDBLIST | xargs -n 1 -P $NPROC sh -c '$AMBER_SCRIPTS_DIR/runAmber_single.sh $1 $NCYC' sh
