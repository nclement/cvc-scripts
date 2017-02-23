#! /bin/bash
# Gets the RMSD between two proteins using TMAlign

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def

#TM_ALIGN=/work/01872/nclement/software/TMalign/TMalign

target=$1
template=$2
$TM_ALIGN $target $template | grep TM-score= | head -n 1 | sed -e 's/.*= //' -e 's/ (if.*//'
