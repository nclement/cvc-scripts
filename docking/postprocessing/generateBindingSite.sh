#! /bin/bash

POST_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR=$POST_SCRIPTS_DIR/../../
source $SCRIPTS_DIR/Makefile.def


########################################################################
# Uses an output from F2Dock and produces a colored image of potential #
# docking locations.                                                   #
########################################################################

inp=$1     # input, should be output from F2Dock, something like dock_0.txt
outp=$2    # output file name, should be .rawc
rmsd=$3    # What RMSD should be considered "contact"? Try something like 5
numConf=$4 # how many configurations should be considered?

if ! grep -q "min RMSD" $inp ; then
  echo "Error: must supply a successful F2Dock output file. Given: [$inp]"
  exit -1
fi

#   staticMoleculeF2d = 1ACB_r_u_samp_6.pdb.scfix_ambermin_1.7.f2d
#   movingMoleculeF2d = 1ACB_l_u_samp_14.pdb.scfix_ambermin.f2d

rec=$( grep -e "staticMoleculeF2d" $inp | sed 's/.* = \(.*\)_1.7.f2d/\1/' )
lig=$( grep -e "movingMoleculeF2d" $inp | sed 's/.* = \(.*\).f2d/\1/' )

# First, generate xforms file.
$POST_SCRIPTS_DIR/extractXForms.pl $numConf $inp > $inp.xforms.txt
$MolSurf_DIR/testBindingSiteUnderPerturbation_F2Dock $rec.rawn $lig.rawn 2 $rmsd $inp.xforms.txt $numConf $outp
