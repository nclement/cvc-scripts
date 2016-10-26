#! /bin/bash
# Take the ZDock output files and create a single csv file with RMSD and ZDOCK Score

# The receptor is needed so we can pull the lines off, and as such should be the
# .zdock.pdb file.
RECEPTOR=$1
LIGAND_ACTUAL=$2
LIG_PATTERN=$3
NUM=$4

DOCKING_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR=${DOCKING_SCRIPTS_DIR}/../
source ${SCRIPTS_DIR}/Makefile.def

CREATE_LIG=${ZDOCK_DIR}/create.pl
#RMSD=${SCRIPTS_DIR}/getRMSD_norot_pymol_zdock.sh
RMSD=${SCRIPTS_DIR}/getRMSD_norot_pymol.sh
  
N_REC_LINES=`wc -l $RECEPTOR | cut -f 1 -d ' '`

echo "rank score rmsd"
for ligand in $LIG_PATTERN; do
  >&2 echo $ligand

  r=$RANDOM
  while [ -f ".$r.energy" ]; do
    r=$RANDOM
  done
  tmpscore=.$r.score
  tmprmsd=.$r.rmsd

  # Generate the conformations and calculate the rmsd
  #echo $CREATE_LIG $ligand $NUM
  $CREATE_LIG $ligand $NUM
  for i in `seq 1 $NUM`; do
    sed -i -e "1,${N_REC_LINES}d;" complex.$i.pdb
    $RMSD $LIGAND_ACTUAL complex.$i.pdb
    \rm complex.$i.pdb
  done > $tmprmsd

  # Get the energy statistics
  numwheader=$(($NUM + 5))
  head -n $numwheader $ligand | tail -n $NUM | tr -s " " | cut -f 7 > $tmpscore
  nl $tmpscore > $tmpscore.tmp

  paste -d ' ' $tmpscore.tmp $tmprmsd

  rm $tmpscore $tmpscore.tmp $tmprmsd
  #echo "press [Enter] to continue"
  #read
done
