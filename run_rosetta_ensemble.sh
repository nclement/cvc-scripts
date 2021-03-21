#! /bin/bash -x
# Rosetta is a bit tricky, but not as bad as F2Dock. This script also helps
# with repeatability going forward. 
#
# Input is a file with ligand conformations, receptor conformations, and output file.
# Won't overrite any of the input files, so it's safe to pass a bunch of files
# in their original form.
#
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
#  cmds.sh 3AAD_f_?_u.fn.txt 3AAD.out
PDBF_l_in=$1
PDBF_r_in=$2
out=$3
#NPROC=${NPROC:-96}
NPROC=96
NUM_DOCK=1000

PDBF_r=${out}.r_fn.txt
PDBF_l=${out}.l_fn.txt

echo "PDB r: $PDBF_r_in, $PDBF_r"
echo "PDB l: $PDBF_l_in, $PDBF_l"
echo "out: $out"

zero_time=$(date +%s)
rm -rf inp && mkdir -p inp
rm -rf output_files && mkdir -p output_files
for fl in $(cat ${PDBF_r_in}); do
  fl_new=$(basename $fl)
  $SCRIPTS_DIR/rosetta_prep_pdb_single.sh $fl inp/${out}.r.${fl_new} R
done
ls inp/${out}.r.*.pdb > ${PDBF_r}
for fl in $(cat $PDBF_l_in); do
  fl_new=$(basename $fl)
  $SCRIPTS_DIR/rosetta_prep_pdb_single.sh $fl inp/${out}.l.${fl_new} L
done
ls inp/${out}.l.*.pdb > ${PDBF_l}

# Grab the first one.
cat $(head -n 1 ${PDBF_r}) $(head -n 1 ${PDBF_l}) > inp/${out}.tog.pdb

one_time=$(date +%s)
echo "############################################"
echo "# RE: Time for initializing is `date -u -d @$(($one_time - $zero_time)) +"%T"`# "
echo "############################################"
# Do some input minimization things. Not parallel, so may as well just run it in one thread.
# This will overwrite the files in ${PFBF_{rl}}, so we can use them again below.
$TACC_ROSETTA_BIN/docking_prepack_protocol.cxx11mpi.linuxiccrelease -database $ROSETTA_DB \
  -ensemble1 ${PDBF_r} -ensemble2 ${PDBF_l} \
  -partners R_L -in:file:s inp/$out.tog.pdb \
  -ex1 -ex2aro \
  -out:path:all output_files -out:suffix _ensemble_prepack -out:overwrite

two_time=$(date +%s)
echo "############################################"
echo "# RE: Time for minimizing inputs is `date -u -d @$(($two_time - $one_time)) +"%T"`# "
echo "############################################"

#ibrun -np $NPROC $TACC_ROSETTA_BIN/docking_protocol.mpiomp.linuxiccrelease -database $ROSETTA_DB -s $tog -dock_pert 3 8 -spin -no_filters -out:overwrite -out:file:scorefile 3AAD.fasc -out:file:fullatom -mute core.io.database -nstruct 1000

# To do global docking, we add the following three options to the options already
# present in global docking.
#
#  -spin
#  -randomize1
#  -randomize2

# Due to the large space sampled, global docking requires a large number of runs
# to converge on a structure, typically 10,000-100,000
NUM_STRUCT=20000
ibrun -np $NPROC $TACC_ROSETTA_BIN/docking_protocol.cxx11mpi.linuxiccrelease -database $ROSETTA_DB \
  -ensemble1 ${PDBF_r} -ensemble2 ${PDBF_l} \
  -dock_pert 3 8 -spin -randomize1 -randomize2 -no_filters \
  -in:file:s inp/$out.tog.pdb \
  -mute core.io.database -ex1 -ex2aro \
  -out:path:all output_files -out:overwrite -out:file:scorefile $out.fasc -out:file:fullatom \
  -nstruct $NUM_STRUCT

three_time=$(date +%s)
echo "############################################"
echo "# RE: Time for doing docking is `date -u -d @$(($three_time - $two_time)) +"%T"`# "
echo "############################################"
# Then create a file with the poses that are in the top NUM_DOCK.
grep -v -e total_score -e "SEQUENCE" output_files/${out}.fasc | sort -k 2n | head -n ${NUM_DOCK} | awk '{print $33}' | sed -e 's~^~output_files/~' -e 's/$/.pdb/' > $out.pdb_fn.txt
