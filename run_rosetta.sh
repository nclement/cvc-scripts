#! /bin/bash
#module use ~/cvc-modules
#module restore rosetta
module load intel/18.0.0  impi/18.0.0
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

$SCRIPTS_DIR/rosetta_prep_pdb.sh $PDB_r $PDB_l $tog

#ibrun -np $NPROC $TACC_ROSETTA_BIN/docking_protocol.mpiomp.linuxiccrelease -database $ROSETTA_DB -s $tog -dock_pert 3 8 -spin -no_filters -out:overwrite -out:file:scorefile 3AAD.fasc -out:file:fullatom -mute core.io.database -nstruct 1000
ibrun -np $NPROC $TACC_ROSETTA_BIN/docking_protocol.cxx11mpi.linuxiccrelease -database $ROSETTA_DB -s $tog -dock_pert 3 8 -spin -no_filters -out:overwrite -out:file:scorefile $out.fasc -out:file:fullatom -mute core.io.database -nstruct 1000
