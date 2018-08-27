#! /bin/bash
#module use ~/cvc-modules
#module restore rosetta
module load rosetta/3.9
module list

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
ADD_TER=$SCRIPTS_DIR/pdb_add_ter.pl
ROSETTA_DB=/work/01872/nclement/software/rosetta_database/3.9

# Can do something like 
#  cmds.sh 3AAD_?_u.pdb 3AAD.out
PDB_l=$1
PDB_r=$2
out=$3
NPROC=${NPROC:-48}

echo "PDB r: $PDB_r"
echo "PDB l: $PDB_l"
echo "out: $out"

tog=$out.pdb
echo " -- tog: $tog"

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

#ibrun -np $NPROC $TACC_ROSETTA_BIN/docking_protocol.mpiomp.linuxiccrelease -database $ROSETTA_DB -s $tog -dock_pert 3 8 -spin -no_filters -out:overwrite -out:file:scorefile 3AAD.fasc -out:file:fullatom -mute core.io.database -nstruct 1000
ibrun -np $NPROC $TACC_ROSETTA_BIN/docking_protocol.cxx11mpi.linuxiccrelease -database $ROSETTA_DB -s $tog -dock_pert 3 8 -spin -no_filters -out:overwrite -out:file:scorefile $out.fasc -out:file:fullatom -mute core.io.database -nstruct 1000
