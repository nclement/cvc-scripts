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
X=${5:-5}
#align_X=10
#align_X=5

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
# Need variables to be defined here.
source $SCRIPTS_DIR/Makefile.def
#PYMOL=/work/01872/nclement/software/pymol/pymol

rms=$($PYMOL -c -p << EOF
load $protR,  pR
load $protRp, pRp
load $protL, pL
load $protLp, pLp
create orig_pL, pLp
alter all, segi=""
alter pR, chain="r"
alter pRp, chain="r"
alter pL, chain="l"
alter pLp, chain="l"
create gold, pL + pR
create testp, pLp + pRp
#select contact_rec_10, ((gold & chain r) & (all within $align_X of (gold and chain l))) & name ca 
#select contact_lig_10, ((gold & chain l) & (all within $align_X of (gold and chain r))) & name ca
select contact_rec, ((gold & chain r) & (all within $X of (gold and chain l))) & name ca 
select contact_lig, ((gold & chain l) & (all within $X of (gold and chain r))) & name ca

# Do the actual alignment of receptors, with no outlier rejection.
#align (testp like contact_rec_10) & name ca, gold & contact_rec_10, cycles=0
align testp & chain r & name ca, gold & contact_rec, cycles=0
# Get the ligand alignments.
align testp & chain l & n. ca, gold & contact_lig, cycles=0, transform=0, object="aln"
#align orig_pL, testp & chain l & aln, object="aln_l"
# Now get the RMSD we care about.
rms_cur testp & aln, gold & aln, matchmaker=-1
#iterate testp & aln, print ID,resi,resn,name,chain
#iterate orig_pL like aln_l, print ID,resi,resn,name,chain
EOF
)
#echo -n "$rms"
if [ $? -eq 0 ]; then
  echo $rms | grep "RMS =" | sed 's/.*=\s\+//' | sed 's/ .*//'
else
  echo "Didn't work. Check the output. Here's what Pymol said:"
  echo -n "$rms"
fi
