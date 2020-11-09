#! /bin/bash

# Will align the CA atoms of one protein with the CA atoms of the other.

prot1=$1
prot2=$2

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#PYMOL=/work/01872/nclement/software/pymol/pymol

rms=$($PYMOL -c -p << EOF
load $prot1, p1
load $prot2, p2
alter all, segi=""
alter all, chain="A"
align name ca and p1, name ca and p2
EOF
)
#echo $rms
echo $rms | grep "RMSD =" | sed 's/.*=\s\+//' | sed 's/ .*//'
