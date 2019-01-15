#! /bin/bash

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
ADD_TER=$SCRIPTS_DIR/pdb_add_ter.pl

PDB_r=$1
PDB_l=$2
tog=$3

# Need to fix the PDB with a few things:
# 1. Get rid of END
# 2. Change residue to be correct.
# 3. Add TER after each chain.
grep "^\(ATOM\|HETATM\|TER\)" $PDB_r | sed 's/^\(.\{21\}\)./\1R/' | $ADD_TER > $tog
grep "^\(ATOM\|HETATM\|TER\)" $PDB_l | sed 's/^\(.\{21\}\)./\1L/' | $ADD_TER >> $tog
# Need to fix it as well.
#$TACC_ROSETTA_TOOLS/protein_tools/scripts/clean_pdb.py $tog ignorechain -database=$ROSETTA_DB
#mv ${tog%.pdb}_ignorechain.pdb $tog
awk -f $SCRIPTS_DIR/rosetta_fixcharmm.awk $tog > ${tog%.pdb}_cleaned.pdb
mv ${tog%.pdb}_cleaned.pdb $tog

