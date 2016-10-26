#! /bin/bash
export SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"

export PDBLIST=$1
export NPROC=$2

cat $PDBLIST | xargs -n 1 -P $NPROC sh -c '$SCRIPTS_DIR/runAmber_single.sh $1' sh
