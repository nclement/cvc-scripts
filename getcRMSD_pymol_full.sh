#! /bin/bash

# Given a protein pair, find the contact-residue RMSD (cRMSD), which is 
# defined as all residues within X\AA of the docked protein complex.
#
# This is done as follows:
#  - Find all residues of A within X\AA of B
#  - Align A with A' along contact RMSD atoms
#  - Compute RMSD over all Ca atoms

protR=$1  # Receptor
protRp=$2 # sampled receptor
protL=$3  # Ligand
protLp=$4 # sampled ligand
# Vreven et al 2015 (zlab5) suggests to use 10\AA as the limit
X=${5:-10}
#align_X=10
#align_X=5

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#PYMOL=/work/01872/nclement/software/pymol/pymol

rms=$($PYMOL -cq $SCRIPTS_DIR/getcRMSD_pymol_full.py -- $protR $protRp $protL $protLp $X)
#echo -n "$rms"
if [ $? -eq 0 ]; then
  #echo $rms
  echo $rms | grep "RMS =" | sed 's/.*=\s\+//'
else
  echo "Didn't work. Check the output. Here's what Pymol said:"
  echo -n "$rms"
fi
