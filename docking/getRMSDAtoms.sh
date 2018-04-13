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
# Vreven et al 2015 (zlab5) suggests to use 10\AA as the limit
distance=${4:-10}


# Some programs to run.
DOCKING_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS=$DOCKING_SCRIPTS_DIR/../
source $DOCKING_SCRIPTS_DIR/../Makefile.def

res=$($PYMOL -c -q -p << EOM
load $RECEPTOR_GOLD, rec_gold
load $LIGAND_GOLD, lig_gold
load $LIGAND_TEST, test_lig

alter all, segi=''

select contact_lig, (lig_gold & (all within $distance of rec_gold) & name ca)
# Print number of atoms.
count_atoms test_lig
count_atoms (test_lig like contact_lig) & name ca
# Need to print atom indices from pAp but locations from pA
my_dict = { 'pA' : [], 'pAp' : [] }
cmd.iterate_state(-1, "(lig_gold & contact_lig)","pA.append((ID,x,y,z,resi,resn,chain))",space=my_dict)
#print my_dict['pA']
cmd.iterate("(test_lig like contact_lig)","pAp.append((ID,resi,resn,chain))",space=my_dict)
#print my_dict['pAp']
for i, x in enumerate(my_dict['pAp']): print "%d %s" % (x[0], "%.3f %.3f %.3f" % tuple(my_dict['pA'][i][1:4]))

# This method doesnt exactly work.
# Print RMSD atoms.
#iterate_state -1, (test_lig like contact_lig) & name ca, print ID,x,y,z
align test_lig like contact_lig, lig_gold & contact_lig, cycles=0, transform=0, object="aln"
rms_cur test_lig & aln, lig_gold & aln, matchmaker=-1
rms_cur test_lig like contact_lig, lig_gold & contact_lig, matchmaker=-1
EOM
)

if [ $? -eq 0 ]; then
  echo "$res"
  printf '%s\n' "$res" | grep -v "^PyMOL" | grep -e " count_atoms" -e "^[^ ]" | sed 's/ count_atoms: \(.*\) atoms/\1/'
else
  echo "Error: Did not work!"
  printf '%s\n' "$res"
fi
