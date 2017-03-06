#! /bin/bash

#############################
# This protocol simply performs F2Dock on a sample of proteins from the inputs.
# 
# Nothing recursive, just specify the number of samples you'd like and the 
# rest will be done.
#############################


RECEPTOR=$1     # The receptor (static--usually bigger of the two)
REC_CHAIN=$2    # Chain identifier(s) of the receptor
REC_GOLD=$3     # The actual receptor (we won't use except for statistics)
LIGAND=$4       # The ligand (moving--usually smaller of the two)
LIG_CHAIN=$5    # Chain identifier(s) of the ligand
LIG_GOLD=$6     # The actual ligand (we won't use except for statistics)
OUTDIR=`readlink -f $7` # Output directory.
NUM_SAMPS=$8    # Number of samples to generate.
stage=${9:-2}

module use ~/cvc-modules
module restore f3dock

# Some other things specified
NPROC=16
NUM_EXTRA=5
NUM_SAMPS_GEND=`echo "$NUM_SAMPS * $NUM_EXTRA" | bc`
RAMACHANDRAN_PROB_FILES="/work/01872/nclement/uq/torsions/Dunbrack/Dunbrack_ctx0_res0.1_fn.txt";
RUNDIR=`pwd`

LIG_SHORT=`basename $LIGAND | sed 's/.pdb//'`
REC_SHORT=`basename $RECEPTOR | sed 's/.pdb//'`

echo "Running F3Dock script with:"
echo " Receptor: $RECEPTOR [$REC_SHORT] chains $REC_CHAIN"
echo " Ligand  : $LIGAND [$LIG_SHORT] chains $LIG_CHAIN"
echo " output directory: $OUTDIR"
echo " number processors: $NPROC"

# Debug tool
#stage=2

SCRIPTS="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS=$SCRIPTS/..
source $SCRIPTS/Makefile.def

DOCK_SCRIPTS=$SCRIPTS/docking

ALIGN_CONTACT=$SCRIPTS/alignContactResidues_pymol.sh
AMBERMIN="$SCRIPTS/amber/runAmber_single.sh"
AMBEREN="$SCRIPTS/amber/getAmberEnergy.sh"
DOCK_BOTH="$DOCK_SCRIPTS/dock_both.sh"
DOCK_BOTH_COARSE="$DOCK_SCRIPTS/dock_both_coarse.sh"
DOCK_LIGAND=$DOCK_SCRIPTS/ligand.sh
#DOCK_SCORE=/work/01872/nclement/scripts/docking/postprocessing/getSortedScores.sh
DOCK_SCORE="$DOCK_SCRIPTS/postprocessing/getSortedPostClusterScores.sh"
FIX_PDB=$SCRIPTS/fix_pdb_residuecount.pl
HINGES_RECURSIVE="$F3DOCK_DIR/FCC_HingeProt"
GET_RMSD=$SCRIPTS/docking/getRMSDAtoms.sh
SAMPLE_RAMACHANDRAN="$F3DOCK_DIR/sampleProtein_Rama"
SINGLES="$SCRIPTS/runSingles_parallel.sh"

echo "############################"
echo "# 0. Generating rmsd atoms #"
echo "############################"
if [ $stage -le 0 ]; then
  mkdir -p $OUTDIR/orig
  # Fix the ligand so it doesn't have any errors.
  $FIX_PDB $LIG_GOLD > $OUTDIR/orig/lig_gold.pdb
  # Then get the .pqr from the ligand.
  cd $OUTDIR/orig
  $DOCK_LIGAND lig_gold.pdb
  cd $RUNDIR

  RMSD_ATOMS=$OUTDIR/orig/rmsd_atoms.txt
  LIG_RMSD=$OUTDIR/orig/lig_gold.pqr
  $PYMOL -qrc $SCRIPTS/get_contact_ca_rmsd_atoms.py -- $LIG_RMSD $REC_GOLD \
        | grep -v "^PyMOL" > $RMSD_ATOMS
fi

echo "###############################################################"
echo "# 1. Performing recursive chain decomposition (One time only) #"
echo "###############################################################"
if [ $stage -le 1 ]; then
  rm -rf $OUTDIR/chains
  mkdir -p $OUTDIR/chains
  cp $LIGAND $OUTDIR/chains
  cp $RECEPTOR $OUTDIR/chains
  cd $OUTDIR/chains
  $HINGES_RECURSIVE $LIGAND   $LIG_CHAIN $OUTDIR/chains F 1 2>&1 | tee $LIGAND.FCC
  $HINGES_RECURSIVE $RECEPTOR $REC_CHAIN $OUTDIR/chains F 1 2>&1 | tee $RECEPTOR.FCC
  cd $RUNDIR
fi

echo "#######################################################"
echo "# 2. Sampling according to Ramachandran distributions #"
echo "#######################################################"
samp_dir=$OUTDIR/samples
rama_args="-R $RAMACHANDRAN_PROB_FILES -N $NUM_SAMPS_GEND --max-clash=10 --max-severe=10 --clash-frac 0.35 --severe-frac 0.35 -s 5"
if [ $stage -le 2 ]; then
  rm -rf $samp_dir
  mkdir -p $samp_dir
  $SAMPLE_RAMACHANDRAN -i $LIGAND   -o $samp_dir/${LIG_SHORT}_samp -S $OUTDIR/chains/${LIG_SHORT}_fluct_resi.txt $rama_args
  $SAMPLE_RAMACHANDRAN -i $RECEPTOR -o $samp_dir/${REC_SHORT}_samp -S $OUTDIR/chains/${REC_SHORT}_fluct_resi.txt $rama_args
fi

echo "#####################################"
echo "# 3. Adding side chains to proteins #"
echo "#####################################"
if [ $stage -le 3 ]; then
  export SCWRL4
  ls $samp_dir/${LIG_SHORT}_samp* $samp_dir/${REC_SHORT}_samp* | xargs -n 1 -P $NPROC sh -c '$SCWRL4 -i $1 -o $1.scfix.pdb' sh
fi

echo "###################################################################"
echo "# 4. Performing (brief) energy minimization and computing energy. #"
echo "###################################################################"
# Go into that directory (AMBER scripts need this to happen)
cd $samp_dir

# Need to save these for future steps of the pipeline.
ls ${LIG_SHORT}_samp*.scfix.pdb > ligands_premin.txt
ls ${REC_SHORT}_samp*.scfix.pdb > recepts_premin.txt

if [ $stage -le 4 ]; then
  # Get each of the files in here and minimize (briefly)
  for file in `cat ligands_premin.txt recepts_premin.txt`; do
    # Only run for 500 iterations.
    $AMBERMIN $file 500
  done
fi

for file in `cat ligands_premin.txt`; do
  echo $file $($AMBEREN $file)
done > ligands_en.txt
for file in `cat recepts_premin.txt`; do
  echo $file $($AMBEREN $file)
done > recepts_en.txt
sort -k 2g ligands_en.txt -o ligands_en.txt
sort -k 2g recepts_en.txt -o recepts_en.txt

# Now need to sort them by energy, then only use the top k, but randomly sort them.
head ligands_en.txt -n $NUM_SAMPS | sort -R > ligands_use.txt
head recepts_en.txt -n $NUM_SAMPS | sort -R > recepts_use.txt


echo "###################################################"
echo "# 5. Aligning proteins back to original receptor. #"
echo "###################################################"
cd $RUNDIR
# Only need to align proteins we'll actually use.
IFS=$'\n' LIG_LIST=($(cat $samp_dir/ligands_use.txt | cut -f 1 -d ' ' | sed 's/.pdb$//'))
IFS=$'\n' REC_LIST=($(cat $samp_dir/recepts_use.txt | cut -f 1 -d ' ' | sed 's/.pdb$//'))
if [ $stage -le 5 ]; then
  # Create these if needs be.
  rm -rf $OUTDIR/aligned
  mkdir $OUTDIR/aligned
  # Copy original here.
  $FIX_PDB $RECEPTOR > $OUTDIR/aligned/orig_receptor.pdb

  for i in `seq 0 $(($NUM_SAMPS - 1))`; do
    lig=${LIG_LIST[$i]}_ambermin.pdb
    rec=${REC_LIST[$i]}_ambermin.pdb
    $FIX_PDB $samp_dir/AMBER/$rec > $OUTDIR/aligned/test_rec.pdb
    $FIX_PDB $samp_dir/AMBER/$lig > $OUTDIR/aligned/$lig # Don't need to change this one.
    $ALIGN_CONTACT $REC_GOLD $LIG_GOLD $OUTDIR/aligned/test_rec.pdb $OUTDIR/aligned/$rec
  done
fi

echo "##########################################"
echo "# 6. Performing F2Dock on coarse samples #"
echo "##########################################"
if [ $stage -le 6 ]; then
  # Might need to clean it.
  if [ -d $OUTDIR/coarse ]; then
    rm -rf $OUTDIR/coarse;
  fi
  mkdir -p $OUTDIR/coarse
fi
cd $OUTDIR/coarse
cp $OUTDIR/samples/ligands_use.txt .
cp $OUTDIR/samples/recepts_use.txt .

# Read the files into an array
IFS=$'\n' LIG_LIST=($(cat ligands_use.txt | cut -f 1 -d ' ' | sed 's/.pdb$//'))
IFS=$'\n' REC_LIST=($(cat recepts_use.txt | cut -f 1 -d ' ' | sed 's/.pdb$//'))
if [ $stage -le 6 ]; then
  # Dock each of the samples.
  for i in `seq 0 $(($NUM_SAMPS - 1))`; do
    lig=${LIG_LIST[$i]}_ambermin.pdb
    rec=${REC_LIST[$i]}_ambermin.pdb
    cp $OUTDIR/aligned/$lig .
    cp $OUTDIR/aligned/$rec .
    # Need to do this to get the .pqr file.
    $DOCK_LIGAND $lig 128
    ( cd $RUNDIR && $GET_RMSD $REC_GOLD $OUTDIR/orig/lig_gold.pqr $OUTDIR/coarse/${lig%.pdb}.pqr > $OUTDIR/coarse/dock_${i}_rmsd_atoms.txt )
    # Include the rmsd atoms just for testing purposes.
    USE_PREV=1 $DOCK_BOTH_COARSE $lig $rec dock_$i.txt dock_${i}_rmsd_atoms.txt
  done
fi

cat dock_*_score.txt | sort -k 2gr > dock_all_score.txt

echo "#######################################"
echo "# 7 Re-Docking will full granulatiry. #"
echo "#######################################"
