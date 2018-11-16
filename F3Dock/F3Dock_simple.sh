#! /bin/bash -x

#############################
# This protocol simply performs F2Dock on a sample of proteins from the inputs.
# 
# Nothing recursive, just specify the number of samples you'd like and the 
# rest will be done.
#############################


TYPE=$1         # The type of docking to perform
RECEPTOR=$2     # The receptor (static--usually bigger of the two)
REC_CHAIN=$3    # Chain identifier(s) of the receptor
REC_GOLD=$4     # The actual receptor (we won't use except for statistics)
LIGAND=$5       # The ligand (moving--usually smaller of the two)
LIG_CHAIN=$6    # Chain identifier(s) of the ligand
LIG_GOLD=$7     # The actual ligand (we won't use except for statistics)
OUTDIR=`readlink -f $8` # Output directory.
NUM_SAMPS=$9    # Number of samples to generate.
NPROC=${10:-16}   # Number of processors to use.
stage=${11:-2}  # Stage to start on.

module use ~/cvc-modules
module restore f3dock

# Some other constants specified.
MAX_HP_LEVEL=4 # Max HingeProt level.
NUM_EXTRA=5
NUM_SAMPS_GEND=`echo "$NUM_SAMPS * $NUM_EXTRA" | bc`
#RAMACHANDRAN_PROB_FILES="/work/01872/nclement/uq/torsions/Dunbrack/Dunbrack_ctx0_res0.1_fn.txt";
RAMACHANDRAN_PROB_FILES="/work/01872/nclement/uq/torsions/Dunbrack/Dunbrack_ctx0_res1_fn.txt";
rama_args="-R $RAMACHANDRAN_PROB_FILES -N $NUM_SAMPS_GEND --max-clash=10 --max-severe=10 --clash-frac 0.35 --severe-frac 0.35 -s 5 --multiply-gauss"
CONTACT_RMSD=5


LIG_SHORT=`basename $LIGAND | sed 's/.pdb//'`
REC_SHORT=`basename $RECEPTOR | sed 's/.pdb//'`

echo "Running F3Dock script with:"
echo " Receptor: $RECEPTOR [$REC_SHORT] chains $REC_CHAIN"
echo " Ligand  : $LIGAND [$LIG_SHORT] chains $LIG_CHAIN"
echo " output directory: $OUTDIR"
echo " number processors: $NPROC"
echo " number of samples: $NUM_SAMPS"

# Debug tool
#stage=2
# stages are:
#   0: Generate RMSD atoms
#   1: Recursive chain decomposition
#   2: Ramachandran sampling
#   3: Adding side chains
#   4: Brief energy minimization
#   5: Align proteins to receptor
#   6: F2Dock (coarse)
#   7: F2Dock (fine)

SCRIPTS="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS=$SCRIPTS/..
source $SCRIPTS/Makefile.def

DOCK_SCRIPTS=$SCRIPTS/docking

ALIGN_CONTACT=$SCRIPTS/alignContactResidues_pymol.sh
AMBERMIN_PAR="$SCRIPTS/amber/runAmber_parallel.sh"
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

one_time=$(date +%s)
## echo "############################"
## echo "# 0. Generating rmsd atoms #"
## echo "############################"
## if [ $stage -le 0 ]; then
##   mkdir -p $OUTDIR/orig
##   # Fix the ligand so it doesnt have any errors.
##   $FIX_PDB $LIGAND > $OUTDIR/$LIG_SHORT.pdb
##   $FIX_PDB $RECEPTOR > $OUTDIR/$REC_SHORT.pdb
##   $FIX_PDB $LIG_GOLD > $OUTDIR/orig/lig_gold.pdb
##   $FIX_PDB $REC_GOLD > $OUTDIR/orig/rec_gold.pdb
##   # Then get the .pqr from the ligand.
##   cd $OUTDIR/orig
##   $DOCK_LIGAND lig_gold.pdb
## 
##   RMSD_ATOMS=$OUTDIR/orig/rmsd_atoms.txt
##   LIG_RMSD=$OUTDIR/orig/lig_gold.pqr
##   $PYMOL -qrc $SCRIPTS/get_contact_ca_rmsd_atoms.py -- $LIG_RMSD $OUTDIR/orig/rec_gold.pdb \
##         | grep -v "^PyMOL" > $RMSD_ATOMS
## fi
## 
two_time=$(date +%s)
echo "############################################"
echo "# Time for step 1 is `date -u -d @$(($two_time - $one_time)) +"%T"`# "
echo "############################################"
echo
echo
echo "###############################################################"
echo "# 1. Performing recursive chain decomposition (One time only) #"
echo "###############################################################"
if [ $stage -le 1 ]; then
  rm -rf $OUTDIR/chains
  mkdir -p $OUTDIR/chains
  cp $OUTDIR/$LIG_SHORT.pdb $OUTDIR/chains
  cp $OUTDIR/$REC_SHORT.pdb $OUTDIR/chains
  cd $OUTDIR/chains
  $HINGES_RECURSIVE $OUTDIR/$LIG_SHORT.pdb $LIG_CHAIN . F $MAX_HP_LEVEL 1 2>&1 | tee $LIG_SHORT.FCC
  $HINGES_RECURSIVE $OUTDIR/$REC_SHORT.pdb $REC_CHAIN . F $MAX_HP_LEVEL 1 2>&1 | tee $REC_SHORT.FCC
fi

three_time=$(date +%s)
echo "############################################"
echo "# Time for step 2 is $(($three_time - $two_time)) #"
echo "############################################"
echo
echo
echo "#######################################################"
echo "# 2. Sampling according to Ramachandran distributions #"
echo "#######################################################"
samp_dir=$OUTDIR/samples
if [ $stage -le 2 ]; then
  rm -rf $samp_dir
  mkdir -p $samp_dir
  cd $OUTDIR
  $SAMPLE_RAMACHANDRAN -i $OUTDIR/$LIG_SHORT.pdb -o $samp_dir/${LIG_SHORT}_samp -S $OUTDIR/chains/${LIG_SHORT}_fluct_resi.txt $rama_args
  $SAMPLE_RAMACHANDRAN -i $OUTDIR/$REC_SHORT.pdb -o $samp_dir/${REC_SHORT}_samp -S $OUTDIR/chains/${REC_SHORT}_fluct_resi.txt $rama_args
fi

four_time=$(date +%s)
echo "############################################"
echo "# Time for step 3 is $(($four_time - $three_time)) #"
echo "############################################"
echo
echo
echo "#####################################"
echo "# 3. Adding side chains to proteins #"
echo "#####################################"
if [ $stage -le 3 ]; then
  export SCWRL4
  ls $samp_dir/${LIG_SHORT}_samp* $samp_dir/${REC_SHORT}_samp* | xargs -n 1 -P $NPROC sh -c '$SCWRL4 -i $1 -o $1.scfix.pdb' sh
fi

five_time=$(date +%s)
echo "############################################"
echo "# Time for step 3 is $(($five_time - $four_time)) #"
echo "############################################"
echo
echo
echo "###################################################################"
echo "# 4. Performing (brief) energy minimization and computing energy. #"
echo "###################################################################"
# Go into that directory (AMBER scripts need this to happen)
cd $samp_dir

# Need to save these for future steps of the pipeline.
ls ${LIG_SHORT}_samp*.scfix.pdb > ligands_premin.txt
ls ${REC_SHORT}_samp*.scfix.pdb > recepts_premin.txt

if [ $stage -le 4 ]; then
  $AMBERMIN_PAR ligands_premin.txt $NPROC 500 # Just do 1 proc for now.
  $AMBERMIN_PAR recepts_premin.txt $NPROC 500
#    # Get each of the files in here and minimize (briefly)
#    for file in `cat ligands_premin.txt recepts_premin.txt`; do
#      # Only run for 500 iterations.
#      $AMBERMIN $file 500
#    done
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


six_time=$(date +%s)
echo "############################################"
echo "# Time for step 4 is $(($six_time - $five_time)) #"
echo "############################################"
echo
echo
echo "###################################################"
echo "# 5. Aligning proteins back to original receptor. #"
echo "###################################################"
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
    $ALIGN_CONTACT $OUTDIR/orig/rec_gold.pdb $OUTDIR/orig/lig_gold.pdb $OUTDIR/aligned/test_rec.pdb $OUTDIR/aligned/$rec
  done
fi

sev_time=$(date +%s)
echo "############################################"
echo "# Time for step 5 is $(($sev_time - $six_time)) #"
echo "############################################"
echo
echo
echo "#####################################################"
echo "# 6. Performing F2Dock on coarse samples (skipping) #" 
echo "#####################################################"
### if [ $stage -le 6 ]; then
###   # Might need to clean it.
###   if [ -d $OUTDIR/coarse ]; then
###     rm -rf $OUTDIR/coarse;
###   fi
###   mkdir -p $OUTDIR/coarse
### fi
### cd $OUTDIR/coarse
### cp $OUTDIR/samples/ligands_use.txt .
### cp $OUTDIR/samples/recepts_use.txt .
### 
### # Read the files into an array
### IFS=$'\n' LIG_LIST=($(cat ligands_use.txt | cut -f 1 -d ' ' | sed 's/.pdb$//'))
### IFS=$'\n' REC_LIST=($(cat recepts_use.txt | cut -f 1 -d ' ' | sed 's/.pdb$//'))
### if [ $stage -le 6 ]; then
###   # Dock each of the samples.
###   for i in `seq 0 $(($NUM_SAMPS - 1))`; do
###     lig=${LIG_LIST[$i]}_ambermin.pdb
###     rec=${REC_LIST[$i]}_ambermin.pdb
###     cp $OUTDIR/aligned/$lig .
###     cp $OUTDIR/aligned/$rec .
###     # Need to do this to get the .pqr file.
###     $DOCK_LIGAND $lig 128
###     $GET_RMSD $OUTDIR/orig/rec_gold.pdb $OUTDIR/orig/lig_gold.pqr $OUTDIR/coarse/${lig%.pdb}.pqr > $OUTDIR/coarse/dock_${i}_rmsd_atoms.txt
### 
### ##    # Include the rmsd atoms just for testing purposes.
### ##    USE_PREV=1 $DOCK_BOTH_COARSE $lig $rec dock_$i.txt dock_${i}_rmsd_atoms.txt
### ##    $DOCK_SCORE dock_$i.txt | sed "s/^/$i /" > dock_${i}_score.txt
###   done
### fi

## # Should be in coarse directory.
## cat dock_*_score.txt | sort -k 2gr > dock_all_score.txt

eight_time=$(date +%s)
echo "############################################"
echo "# Time for step 6 is $(($eight_time - $sev_time)) #"
echo "############################################"
echo
echo
echo "#######################################"
echo "# 7. Re-Docking will full granulatiry. #"
echo "#######################################"
if [ $stage -le 7 ]; then
  # Might need to clean it.
  if [ -d $OUTDIR/fine ]; then
    rm -rf $OUTDIR/fine;
  fi
  mkdir -p $OUTDIR/fine
fi
cd $OUTDIR/fine
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
    # Need to do this to get the rmsd_atoms file.
    $DOCK_LIGAND $lig 128
    #$GET_RMSD $OUTDIR/orig/rec_gold.pdb $OUTDIR/orig/lig_gold.pqr ${lig%.pdb}.pqr $CONTACT_RMSD > dock_${i}_rmsd_atoms.txt

    # Include the rmsd atoms just for testing purposes.
    #USE_PREV=1 $DOCK_BOTH $lig $rec dock_$i.txt dock_${i}_rmsd_atoms.txt
    #NUM_THREADS=$NPROC TYPE=$TYPE $DOCK_BOTH $lig $rec dock_$i.txt dock_${i}_rmsd_atoms.txt
    NUM_THREADS=$NPROC TYPE=$TYPE $DOCK_BOTH $lig $rec dock_$i.txt
    $DOCK_SCORE dock_$i.txt | sed "s/^/$i /" > dock_${i}_score.txt
  done
fi

# Should be in fine directory.
cat dock_*_score.txt | sort -k 2gr > dock_all_score.txt

e_time=$(date +%s)
echo "############################################"
echo "# Time for step 7 is $(($e_time - $eight_time)) #"
echo "############################################"
echo
echo
