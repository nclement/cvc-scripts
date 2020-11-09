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
NPROC=${10:-16} # Number of processors to use.
stage=${11:-2}  # Stage to start on.

PROTOCOL=${PROTOCOL:-f3d}  # Sampling protocol to use; should be `f3d`, `atom`, or `hinge`.
F2D_PROC=${F2D_PROC:-$NPROC}

module use ~/cvc-modules
module restore f3dock

# Some other constants specified.
MAX_HP_LEVEL=4 # Max HingeProt level.
NUM_EXTRA=5
NUM_SAMPS_GEND=`echo "$NUM_SAMPS * $NUM_EXTRA" | bc`
# Default std of 5 unless specified
R_stdv=${R_stdv:-5}
# Context for Ramachandran distributions.
R_ctx=${R_ctx:-1}
CONTACT_RMSD=5


LIG_SHORT=`basename $LIGAND | sed 's/.pdb//'`
REC_SHORT=`basename $RECEPTOR | sed 's/.pdb//'`
samp_dir=$OUTDIR/samples
# Need to compute these offline.
#hinges_prefix=$OUTDIR/chains/
hinges_prefix=/work/01872/nclement/F2dock_rerun/benchmark5/hinges/hp.5.

echo "Running F3Dock script with:"
echo " Protocol: $PROTOCOL"
echo " Receptor: $RECEPTOR [$REC_SHORT] chains $REC_CHAIN"
echo " Ligand  : $LIGAND [$LIG_SHORT] chains $LIG_CHAIN"
echo " output directory: $OUTDIR"
echo " hinges dir: $hinges_prefix"
echo " number processors: $NPROC"
echo " number of samples: $NUM_SAMPS"

# Debug tool
#stage=2
# stages are:
#   0: Generate RMSD atoms - don't do this anymore
#   1: Recursive chain decomposition - do this offline
#   2: Ramachandran sampling
#   3: Adding side chains
#   4: Brief energy minimization
#   5: Align proteins to receptor
#   6: F2Dock (coarse)
#   7: F2Dock (fine)

F3DOCK_SCRIPTS="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS=$F3DOCK_SCRIPTS/..
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




six_time=$(date +%s)
# These have all been exported to conformation_gen.sh
export R_stdv R_ctx
$F3DOCK_SCRIPTS/conformation_gen.sh $PROTOCOL $OUTDIR $LIGAND $RECEPTOR $NUM_SAMPS $NPROC




sev_time=$(date +%s)
echo "############################################"
echo "# Time for step 5 is `date -u -d @$(($sev_time - $six_time)) +"%T"`# "
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
echo "# Time for step 6 is `date -u -d @$(($eight_time - $sev_time)) +"%T"`# "
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
    ## # Need to do this to get the rmsd_atoms file.
    ## #$DOCK_LIGAND $lig 128
    ## #$GET_RMSD $OUTDIR/orig/rec_gold.pdb $OUTDIR/orig/lig_gold.pqr ${lig%.pdb}.pqr $CONTACT_RMSD > dock_${i}_rmsd_atoms.txt
    ## 
    ## # Include the rmsd atoms just for testing purposes.
    ## #USE_PREV=1 $DOCK_BOTH $lig $rec dock_$i.txt dock_${i}_rmsd_atoms.txt
    ## #NUM_THREADS=$NPROC TYPE=$TYPE $DOCK_BOTH $lig $rec dock_$i.txt dock_${i}_rmsd_atoms.txt
    RES_CONT_FILTER=true HBOND_FILTER=false NUM_THREADS=$F2D_PROC TYPE=$TYPE \
      $DOCK_BOTH $lig $rec dock_$i.txt
    $DOCK_SCORE dock_$i.txt | sed "s/^/$i /" > dock_${i}_score.txt
  done
fi

# Should be in fine directory.
cat dock_*_score.txt | sort -k 2gr > dock_all_score.txt

e_time=$(date +%s)
echo "############################################"
echo "# Time for step 7 is `date -u -d @$(($e_time - $eight_time)) +"%T"`# "
echo "############################################"
echo
echo
