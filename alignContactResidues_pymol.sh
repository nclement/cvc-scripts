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
X=${5:-10} # Vreven et al 2015 (zlab5) suggests to use 10\AA as the limit

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#PYMOL=/work/01872/nclement/software/pymol/pymol

cmd=$($PYMOL -c -p << EOF
load $protR, pR
load $protL, pL
load $protA, pA
alter all, segi=""
# select all CA atoms in receptor that have any atom within $X of any atom in ligand
select contact, (pR & (all within $X of pL)) & name ca
  # selects residues from pA that have names and resi atoms match pR
#my_dict = { 'pR' : [] }
#cmd.iterate("(contact & pR)","pR.append((resi,resn,chain))",space=my_dict)
#print my_dict['pR']
#my_dict = { 'pA' : [] }
#cmd.iterate("(pA like contact)","pA.append((resi,resn,chain))",space=my_dict)
#print my_dict['pA']
#align contact and pR, pA like contact
# Avoids mismatches in residue numbering, etc.
align pA like contact, contact and pR
save $output, pA
EOF
)
#echo $cmd
if [ $? -eq 0 ]; then
  echo "Success!"
else
  echo "Didn't work. Check output. Here's what Pymol said:"
  echo $cmd
fi
