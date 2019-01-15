#! /bin/bash

# Given a protein pair, find the contact-residue RMSD (cRMSD), which is 
# defined as all residues within X\AA of the docked protein complex.

protA=$1  # Original protein A
protB=$2  # Docked protein B
# Vreven et al 2015 (zlab5) suggests to use 10\AA as the limit
X=${3:-10}

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#PYMOL=/work/01872/nclement/software/pymol/pymol

cRMSD=$($PYMOL -c -p << EOF
load $protA,  pA
load $protB, pB
alter all, segi=""
alter pA, chain="A"
select contact, name ca & pA & (all within $X of pB)
my_dict = { 'pA' : [] }
cmd.iterate("(contact & pA)","pA.append(resi)",space=my_dict)
print "contact RMSD: ", my_dict['pA']
EOF
)
#echo $cRMSD
echo $cRMSD | grep -o "contact RMSD: \\[.*" | sed -e 's/.*\[//' -e 's/\].*//' -e "s/'//g" | sed -e "s/, /+/g"
