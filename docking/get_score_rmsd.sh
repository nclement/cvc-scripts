#! /bin/bash
# Calculates the RMSD and F2Dock score for each individual docked ligand

DOCKING_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"

F2DOCK_2_XFORMS=${DOCKING_SCRIPTS_DIR}/postprocessing/converter
CONFORMATION_GEN=${DOCKING_SCRIPTS_DIR}/postprocessing/conformationGenerator
#RMSD=/work/01872/nclement/scripts/getRMSD_TM.sh
RMSD=${DOCKING_SCRIPTS_DIR}/../getRMSD_norot_pymol.sh

DOCKING_FOLDER=$1 # The folder where the docking output files reside
DOCKING_EXT=$2    # The extension of the docking output files
NPOSES=$3         # Number of top poses to report
LIGAND_ACTUAL=$4  # Optional, the ligand to compare the RMSD to

# For the latest version of F2dock, it's this:
RANK_COL=2
SCORE_COL=3
SCALE_COL=31
RMSD_COL=30

echo "rank score scale rmsd"
if [[ "$#" -eq "4" ]]; then
  LIGAND_ACTUAL=$4
  # Need to do some postprocessing
  for file in $DOCKING_FOLDER/*$DOCKING_EXT; do
    >&2 echo $file
    # Create a temp file that we'll delete afterwards
    r=$RANDOM
    while [ -f ".$r.energy" ]; do
      r = $RANDOM
    done
    tmpscore=.$r.score
    tmprmsd=.$r.rmsd

    # Generate the xforms and calculate the RMSD
    # Need one less xform
    nminusone=$(($NPOSES-1))
    $F2DOCK_2_XFORMS $file $file.xforms $file.rmsds &>/dev/null
    $CONFORMATION_GEN ${file/pdb*/pdb} $file.xforms 0 $nminusone &>/dev/null
    for conf in ${file/.pdb*/_conf}*; do
      $RMSD $LIGAND_ACTUAL $conf
    done > $tmprmsd
    # Need to reverse this file
    tac $tmprmsd > $tmprmsd.tmp
    mv $tmprmsd.tmp $tmprmsd
   
    # Generate the energy statistics
    scale=`grep "score scale down" $file | sed 's/.*= //'`
    grep -B $NPOSES "END PEAKS" $file | head -n $NPOSES | tr -s " " | sed "s/$/ $scale/" | cut -f $RANK_COL,$SCORE_COL,$SCALE_COL -d " " > $tmpscore
    paste -d ' ' $tmpscore $tmprmsd

    # Clean up tmp files
    rm $tmpscore $tmprmsd
  done
else
  # All the information is in the docking file
  for file in $DOCKING_FOLDER/*$DOCKING_EXT; do
    scale=`grep "score scale down" $file | sed 's/.*= //'`
    grep -B $NPOSES "END PEAKS" $file | head -n $NPOSES | tr -s " " | sed "s/$/ $scale/" | cut -f $RANK_COL,$SCORE_COL,$SCALE_COL,$RMSD_COL -d " "
  done
fi
