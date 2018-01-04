#! /bin/bash

##########################################################
# Creating an rmsd_atoms files for F2Dock can be kind of #
# tedious because the atom numberings must match in the  #
# gold standard and the target model. This script will   #
# generate it for you, and ensure that it will work with #
# F2Dock.                                                #
#                                                        #
# As input, pass in the (actual) receptor and ligand     #
# files as well as the template ligand to dock with the  #
# receptor.                                              #
#                                                        #
# Something like:                                        #
#   getRMSDAtoms.sh 1ATN_r_b.pdb 1ATN_l_b.pdb 1ATN_l_u.pdb > 1ATN_rmsd_atoms.txt
##########################################################

RECEPTOR_GOLD=$1
LIGAND_GOLD=$2
LIGAND_TEST=$3


# Some programs to run.
DOCKING_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS=$DOCKING_SCRIPTS_DIR/../
source $DOCKING_SCRIPTS_DIR/../Makefile.def

#echo $PYMOL -qrc $SCRIPTS/get_contact_ca_rmsd_atoms.py -- $LIGAND_GOLD $RECEPTOR_GOLD

count=1
IFS=$'\n'
$PYMOL -qrc $SCRIPTS/get_contact_ca_rmsd_atoms.py -- $LIGAND_GOLD $RECEPTOR_GOLD \
  | grep -v "^PyMOL" | \
  ( read numlines
    echo $(grep -c "^ATOM" $LIGAND_TEST)
    read numthings
    echo $numthings
  while read line; do
    atomnum=`echo $line | sed 's/ .*//'`
    location=`echo $line | sed 's/^[0-9]\+ //'`
    former=$(grep -m 1 "ATOM \\+$atomnum " $LIGAND_GOLD | sed 's/^.\{11\}\(.\{15\}\).*/\1/')
    # This is a terrible hack
    if ! grep -q "ATOM.*$former" $LIGAND_TEST; then
      echo >&2 "ERROR: Could not find residue!!! looking for [$former]"
      # Strip off the residue name.
      former=$(echo $former | sed 's/^\(.\{6\}\).../\1.../')
    fi
    current=$(grep -m 1 "ATOM.*$former" $LIGAND_TEST | sed 's/^ATOM \+\([^ ]\+\).*/\1/')
    echo >&2 $atomnum $former becomes [$current] $location
    echo $current $location
  done )
