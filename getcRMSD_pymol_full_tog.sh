#! /bin/bash

# Given a protein pair, find the contact-residue RMSD (cRMSD), which is 
# defined as all residues within X\AA of the docked protein complex.
#
# This is done as follows:
#  - Find all residues of A within X\AA of B
#  - Align A with A' along contact RMSD atoms
#  - Compute RMSD over all Ca atoms
#
# This script accepts two proteins that have both chains in them and computes
# the cRMSD of ligand-ligand only.

prot=$1  # Protein
samp=$2  # Sampled protein
chainsL=$3
chainsR=$4
X=${5:-5}
#align_X=10
#align_X=5

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#PYMOL=/work/01872/nclement/software/pymol/pymol

rms=$($PYMOL -qc -p << EOF
print('chains are R:[$chainsR] and L:[$chainsL]')
load $prot,  gold
load $samp, testp
alter all, segi=""
select contact_rec, ((gold & chain $chainsR) & (all within $X of (gold and chain $chainsL))) & name ca 
select contact_lig, ((gold & chain $chainsL) & (all within $X of (gold and chain $chainsR))) & name ca
select contact_both, contact_rec + contact_lig
# Do the actual alignment with no outlier rejection.
align testp, gold & contact_all, cycles=0
EOF
)
#echo -n "$rms"
if [ $? -eq 0 ]; then
  echo $rms | grep "RMS =" | sed 's/.*=\s\+//' | sed 's/ .*//'
  #echo $rms | grep "RMS =" | sed 's/.*=\s\+//'
else
  echo "Didn't work. Check the output. Here's what Pymol said:"
  echo -n "$rms"
fi
