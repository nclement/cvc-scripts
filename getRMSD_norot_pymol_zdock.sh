#! /bin/bash
# ZDOCK does some weird things to these proteins, so we need to handle
# the complex (2nd arg) differently than the normal PDB.
prot1=$1
complex=$2


SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#PYMOL=/work/01872/nclement/software/pymol/pymol

if [[ "$#" -eq "3" ]]; then
  echo "$PYMOL -c -p
  load $prot1, p1
  load $complex, p2
  remove p2 and chain A
  alter p1, chain="A"
  alter p2, chain="A"
  renumber p1
  renumber p2
  rms_cur name ca and p1, name ca and p2"
fi

rms=$($PYMOL -c -p << EOF
load $prot1, p1
load $complex, p2
remove p2 and chain A
alter p1, chain="A"
alter p2, chain="A"
renumber p1
renumber p2
rms_cur name ca and p1, name ca and p2
EOF
)
echo $rms | grep "RMS =" | sed 's/.*=\s\+//' | sed 's/ .*//'
