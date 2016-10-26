#! /bin/bash

prot1=$1
prot2=$2

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#PYMOL=/work/01872/nclement/software/pymol/pymol
rms=$($PYMOL -c -p << EOF
load $prot1, p1
load $prot2, p2
alter p1, chain="A"
alter p2, chain="A"
rms_cur name ca and p1, name ca and p2
EOF
)
#echo $rms
echo $rms | grep "RMS =" | sed 's/.*=\s\+//' | sed 's/ .*//'
