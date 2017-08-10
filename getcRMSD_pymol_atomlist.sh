#! /bin/bash

# Given a single protein and a list of contact residues, compute the RMSD
# between the protein and an additional sampled protein at only the contact
# residues.
#
# Uses output from getcRMSD_atoms.sh.
#
protA=$1  # Original protein
protAp=$2 # Sampled protein
contact=$3  # Contact residues (from e.g. getcRMSD_atoms.sh)

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#PYMOL=/work/01872/nclement/software/pymol/pymol

rms=$($PYMOL -c -p << EOF
load $protA,  pA
load $protAp, pAp
alter all, segi=""
alter pA, chain="A"
alter pAp, chain="A"
select contact, name ca & pA & resi $contact
#my_dict = { 'pA' : [] }
#cmd.iterate("(contact & pA)","pA.append((resi,resn,chain))",space=my_dict)
#print my_dict['pA']
#my_dict = { 'pAp' : [] }
#select thing, pAp l. contact
#cmd.iterate("(pAp like contact)","pAp.append((resi,resn,chain))",space=my_dict)
#print my_dict['pAp']
# rms doesnt do the actual fitting, just finds the rmsd.
rms contact and pA, pAp like contact
EOF
)
#echo $rms
#echo $rms | grep "RMS =" | sed 's/.*=\s\+//' | sed 's/ atoms.*/ atoms)/'
echo $rms | grep "RMS =" | sed 's/.*=\s\+//' | sed 's/ .*//'
