#! /bin/bash

# Given a protein pair, find the contact-residue RMSD (cRMSD), which is 
# defined as all residues within X\AA of the docked protein complex.
#
# This is done as follows:
#  - Find all residues of A within X\AA of B
#  - Align A with A'
#  - Compute RMSD over all Ca atoms

protA=$1  # Original protein
protAp=$2 # Sampled protein
protB=$3  # Docked protein
X=${4:-5}

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#PYMOL=/work/01872/nclement/software/pymol/pymol

rms=$($PYMOL -c -p << EOF
load $protA,  pA
load $protAp, pAp
load $protB, pB
alter all, segi=""
alter pA, chain="A"
alter pAp, chain="A"
select contact, name ca & pA & (all within $X of pB)
#my_dict = { 'pA' : [] }
#cmd.iterate("(contact & pA)","pA.append((resi,resn,chain))",space=my_dict)
#print my_dict['pA']
#my_dict = { 'pAp' : [] }
#select thing, pAp l. contact
#cmd.iterate("(pAp like contact)","pAp.append((resi,resn,chain))",space=my_dict)
#print my_dict['pAp']
align contact and pA, pAp like contact
EOF
)
#echo $rms
echo $rms | grep "RMS =" | sed 's/.*=\s\+//' | sed 's/ .*//'
