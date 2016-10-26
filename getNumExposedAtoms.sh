#! /bin/bash

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#MolSurf=/work/01872/nclement/MolSurf/release/bin/MolSurf

PDB=$1

tmpfile=/tmp/surfaceAtoms.txt
$MolSurf -surfaceAtoms $PDB $tmpfile &> /dev/null
grep -c "^E" $tmpfile
