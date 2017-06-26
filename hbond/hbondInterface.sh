#! /bin/bash

HBOND_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR="$HBOND_SCRIPTS_DIR/../"
source $SCRIPTS_DIR/Makefile.def

EXTRACT_INTERFACE=$MOLSURF_BIN/testExtractInterfaceAtoms
GET_NUM_INTERACT=$MOLSURF_BIN/getNumAtomsOnInterface

# Computes all the atoms on the interface, and reports their atomic energies.
PDBA=$1
PDBB=$2
OUT=$3

# Get the chains for each. Charmm ignores HETATMs, so we can, too
CHAINS_A=$(cat $PDBA | grep -e "^ATOM" | cut -c 22 | uniq | tr -d '\n')
CHAINS_B=$(cat $PDBB | grep -e "^ATOM" | cut -c 22 | uniq | tr -d '\n')

# Get the final residue numbers
N_RESI_A=$(grep -e "^ATOM" $PDBA | cut -c 23-26 | tail -n 1)
N_RESI_B=$(grep -e "^ATOM" $PDBB | cut -c 23-26 | tail -n 1)

# Concatenate them together.
OUT_TOG=$OUT-tog.pdb
$SCRIPTS_DIR/fix_pdb_atomcount.pl <(cat $PDBA $PDBB) > $OUT_TOG
sed -i "1iREMARK Generated by a script for hbondInterface" $OUT_TOG
sed -i "1iHBOND $OUT-tog.hbond" $OUT_TOG

# Get the hbonded atoms
echo $HBOND_SCRIPTS_DIR/runHbondSingle.sh $OUT_TOG
$HBOND_SCRIPTS_DIR/runHbondSingle.sh $OUT_TOG

OUT_MIN=$OUT-tog_cmin.pdb

# Then, get the correct interface atoms
# First, split the minimized structure into chains
$SCRIPTS_DIR/extract_chains.pl $OUT_MIN
# Then, concatenate them together to the respective files.
cat ${OUT_MIN%.pdb}-[$CHAINS_A].pdb > $OUT-tog_cmin_chainsA.pdb
cat ${OUT_MIN%.pdb}-[$CHAINS_B].pdb > $OUT-tog_cmin_chainsB.pdb
# Then, get the interface atoms
echo $EXTRACT_INTERFACE $OUT-tog_cmin_chainsA.pdb $OUT-tog_cmin_chainsB.pdb ${PDBA%.pdb}-int.txt ${PDBB%.pdb}-int.txt
$EXTRACT_INTERFACE $OUT-tog_cmin_chainsA.pdb $OUT-tog_cmin_chainsB.pdb ${PDBA%.pdb}-int.txt ${PDBB%.pdb}-int.txt
cat ${PDBA%.pdb}-int.txt ${PDBB%.pdb}-int.txt > ${OUT}_int.txt

# Get the number of atoms in either.
N_A=$(grep -e "^ATOM" $OUT_MIN | cut -c 22 | grep -c "[$CHAINS_A]")
N_B=$(grep -e "^ATOM" $OUT_MIN | cut -c 22 | grep -c "[$CHAINS_B]")
echo "Number of atoms in A and B is [$N_A] [$N_B] from chains [$CHAINS_A] [$CHAINS_B]"

# Output is in a file called $OUT-tog.hbond. Need to add a single line at
# the beginning that tells how many atoms are in the receptor and ligand.
sed -i "1iREMARK Generated by a script for hbondInterface" $OUT_MIN
# Get the last REMARK line for the HBOND line
ln=$(( $(grep -n REMARK $OUT_MIN | tail -n 1 | cut -f 1 -d ':') + 1))
sed -i "${ln}iHBOND $OUT.hbond" $OUT_MIN
sed -i "1i$N_A $N_B" $OUT-tog.hbond
cp $OUT_MIN $OUT.pdb
cp $OUT-tog.hbond  $OUT.hbond
cp $OUT-tog_cmin.mol2 $OUT.mol2

# Get the number of atoms on the interface.
echo Interface counts 2AA: $($GET_NUM_INTERACT $OUT-tog_cmin_chainsA.pdb $OUT-tog_cmin_chainsB.pdb 2 0)
echo Interface counts 5AA: $($GET_NUM_INTERACT $OUT-tog_cmin_chainsA.pdb $OUT-tog_cmin_chainsB.pdb 5 0)

# Will create three files. First two are like a PQR file, as follows
#
#   ATOM_NUM_1 E_COUL E_LJ E_HB
#   ATOM_NUM_2 E_COUL E_LJ E_HB
#   ...
#   ATOM_NUM_N E_COUL E_LJ E_HB
#
# The same file will be created for each protein.
# 
# Then an additional interaction file will be created, as follows:
# ATOM_NUM_A ATOM_NUM_B RES_NUM_

## if [ -z ${NO_CLEANUP-x} ]; then
##   # Cleanup
##   rm $OUT_TOG
##   rm ${PDBA%.pdb}-int.txt ${PDBB%.pdb}-int.txt
##   rm $OUT-tog*
## else
##   echo "Not cleaning up (because NO_CLEANUP set)"
## fi
## 
