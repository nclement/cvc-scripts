#! /bin/bash
# Rosetta is a bit tricky, but not as bad as F2Dock. This script also helps
# with repeatability going forward. 
#
# Input is a ligand, receptor, and output file.
# Output is three items:
#  1) a file ${out}.fasc: the scores of each decoy (see Rosetta manual for more details)
#  2) a file ${out}.pdb_fn.txt: the names of the PDBs of the top ${NUM_DOCK} poses
#     -> will dock 20000 different decoys, then this file only has 1000
#  3) all the 20k different poses, in the form of ${out}_XXXX.pdb
#
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
NUM_DOCK=1000

echo "PDB r: $PDB_r"
echo "PDB l: $PDB_l"
echo "out: $out"

tog=$out.pdb
echo " -- tog: $tog"

$SCRIPTS_DIR/rosetta_prep_pdb.sh $PDB_r $PDB_l $tog

#ibrun -np $NPROC $TACC_ROSETTA_BIN/docking_protocol.mpiomp.linuxiccrelease -database $ROSETTA_DB -s $tog -dock_pert 3 8 -spin -no_filters -out:overwrite -out:file:scorefile 3AAD.fasc -out:file:fullatom -mute core.io.database -nstruct 1000

# To do global docking, we add the following three options to the options already
# present in global docking.
#
#  -spin
#  -randomize1
#  -randomize2

# Due to the large space sampled, global docking requires a large number of runs
# to converge on a structure, typically 10,000-100,000
NUM_STRUCT=10000
ibrun -np $NPROC $TACC_ROSETTA_BIN/docking_protocol.cxx11mpi.linuxiccrelease -database $ROSETTA_DB -s $tog -dock_pert 3 8 -spin -randomize1 -randomize2 -no_filters -out:overwrite -out:file:scorefile $out.fasc -out:file:fullatom -mute core.io.database -ex1 -ex2aro -nstruct $NUM_STRUCT

# Then create a file with the poses that are in the top NUM_DOCK.
grep -v -e total_score -e "SEQUENCE" ${out}.fasc | sort -k 2n | head -n ${NUM_DOCK} | awk '{print $27}' | sed 's/$/.pdb/' > $out.pdb_fn.txt
