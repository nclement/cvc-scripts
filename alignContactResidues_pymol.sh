#! /bin/bash

# Given a protein pair (receptor and ligand, R and L) and a third protein, A,
# align the residues from A to the residues on R that are in contact with L,
# defined as all residues within X\AA of R and L.
#
# This is done as follows:
#  - Find all residues of R within X\AA of L
#  - Align R with L
#  - Output A'
#
# Note that A should have the same residue numbering as R for this to work.

protR=$1  # Receptor
protL=$2  # Ligand
protA=$3  # Protein A
output=$4 # Output protein A'
X=${5:-5}

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#PYMOL=/work/01872/nclement/software/pymol/pymol

rms=$($PYMOL -c -p << EOF
load $protR, pR
load $protL, pL
load $protA, pA
alter all, segi=""
select contact, name ca & pR & (all within $X of pL)
  # selects residues from pA that have names and resi atoms match pR
align contact and pR, pA like contact
save $output, pA
EOF
)
echo $rms
