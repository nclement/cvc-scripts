#! /bin/bash -x

F3DOCK_SCRIPTS="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS=$F3DOCK_SCRIPTS/..
source $SCRIPTS/Makefile.def

FIX_PDB_A=$SCRIPTS/fix_pdb_atomcount.pl
FIX_PDB=$SCRIPTS/fix_pdb_residuecount.pl
SAMPLE_RAMACHANDRAN="$F3DOCK_DIR/sampleProtein_Rama"
AMBERMIN_PAR="$SCRIPTS/amber/runAmber_parallel.sh"
AMBEREN="$SCRIPTS/amber/getAmberEnergy.sh"
# Need to compute these offline.
#hinges_prefix=$OUTDIR/chains/
hinges_prefix=/work/01872/nclement/F2dock_rerun/benchmark5/hinges/hp.5.

PROTOCOL=$1  # Sampling protocol to use; should be `f3d`, `atom`, or `hinge`.
OUTDIR=$2
LIGAND=$3
RECEPTOR=$4
NUM_SAMPS=$5
NPROC=$6


# Here we define some constants
NUM_EXTRA=5
NUM_SAMPS_GEND=`echo "$NUM_SAMPS * $NUM_EXTRA" | bc`
LIG_SHORT=`basename $LIGAND | sed 's/.pdb//'`
REC_SHORT=`basename $RECEPTOR | sed 's/.pdb//'`


samp_dir=$OUTDIR/samples
rm -rf $samp_dir && mkdir -p $samp_dir


gen_atom_samples() {
  # This wrecks havoc with the whole process, so clear all these things out.
  # Also fix atom count (because F2Dock has issues with it), but don't touch the
  # residue numbering.
  sed $LIGAND   -e '/^REMARK 350/d' -e '/^HETATM/d' | $FIX_PDB_A false > $samp_dir/${LIG_SHORT}.cleaned.pdb
  sed $RECEPTOR -e '/^REMARK 350/d' -e '/^HETATM/d' | $FIX_PDB_A false > $samp_dir/${REC_SHORT}.cleaned.pdb
  # Generate samples first.
  if [[ $(wc -l ligands_premin.txt | awk '{print $1}') -lt $NUM_SAMPS_GEND ]]; then
    $SCRIPTS/F3Dock/sample_atomic.sh $samp_dir/${LIG_SHORT}.cleaned.pdb $NUM_SAMPS_GEND $samp_dir/${LIG_SHORT}
  fi
  if [[ $(wc -l recepts_premin.txt | awk '{print $1}') -lt $NUM_SAMPS_GEND ]]; then
    $SCRIPTS/F3Dock/sample_atomic.sh $samp_dir/${REC_SHORT}.cleaned.pdb $NUM_SAMPS_GEND $samp_dir/${REC_SHORT}
  fi

  # Save these for future steps of the pipeline.
  (cd $samp_dir;
   ls ${LIG_SHORT}_samp*.pdb > ligands_premin.txt;
   ls ${REC_SHORT}_samp*.pdb > recepts_premin.txt;)
}

gen_hinge_samples() {
  # This wrecks havoc with the whole process, so clear all these things out.
  sed $LIGAND   -e '/^REMARK 350/d' -e '/^HETATM/d' | $FIX_PDB_A false > $samp_dir/${LIG_SHORT}.cleaned.pdb
  sed $RECEPTOR -e '/^REMARK 350/d' -e '/^HETATM/d' | $FIX_PDB_A false > $samp_dir/${REC_SHORT}.cleaned.pdb

  # First, convert from hinges into the thing we like.
  # Need to generate >>NUM_SAMPS_GEND because we may have clashes. But then filter later if we need to.
  nsamp_over=$(($NUM_SAMPS_GEND*100))
  # Then generate samples.
  if [[ $(wc -l ligands_premin.txt | awk '{print $1}') -lt $NUM_SAMPS_GEND ]]; then
    $F3DOCK_SCRIPTS/HingeProt_to_hinges.pl ${hinges_prefix}${LIG_SHORT}_fcc.txt 1 $nsamp_over > $samp_dir/${LIG_SHORT}_hinges.txt
    $F3DOCK_DIR/sampleHinges $samp_dir/${LIG_SHORT}.cleaned.pdb $samp_dir/${LIG_SHORT}_hinges.txt $samp_dir/${LIG_SHORT} 1 $NPROC
  fi
  if [[ $(wc -l recepts_premin.txt | awk '{print $1}') -lt $NUM_SAMPS_GEND ]]; then
    $F3DOCK_SCRIPTS/HingeProt_to_hinges.pl ${hinges_prefix}${REC_SHORT}_fcc.txt 1 $nsamp_over > $samp_dir/${REC_SHORT}_hinges.txt
    $F3DOCK_DIR/sampleHinges $samp_dir/${REC_SHORT}.cleaned.pdb $samp_dir/${REC_SHORT}_hinges.txt $samp_dir/${REC_SHORT} 1 $NPROC
  fi

  # Save these for future steps of the pipeline.
  (cd $samp_dir;
   ls ${LIG_SHORT}_N*.pdb | head -n $NUM_SAMPS_GEND > ligands_premin.txt;
   ls ${REC_SHORT}_N*.pdb | head -n $NUM_SAMPS_GEND > recepts_premin.txt;)
}

gen_f3d_samples() {
  sed $LIGAND   -e '/^REMARK 350/d' -e '/^HETATM/d' | $FIX_PDB_A false > $samp_dir/${LIG_SHORT}.cleaned.pdb
  sed $RECEPTOR -e '/^REMARK 350/d' -e '/^HETATM/d' | $FIX_PDB_A false > $samp_dir/${REC_SHORT}.cleaned.pdb

  # Default std of 5 unless specified
  R_stdv=${R_stdv:-5}
  # Context for Ramachandran distributions.
  R_ctx=${R_ctx:-1}
  #RAMACHANDRAN_PROB_FILES="/work/01872/nclement/uq/torsions/Dunbrack/Dunbrack_ctx01_res0.1_fn.txt";
  RAMACHANDRAN_PROB_FILES="/work/01872/nclement/uq/torsions/Dunbrack/Dunbrack_ctx01_res1_fn.txt";
  rama_args="-R $RAMACHANDRAN_PROB_FILES -N $NUM_SAMPS_GEND -C ${R_ctx} --max-clash=10 --max-severe=10 --clash-frac 0.35 --severe-frac 0.35 -s $R_stdv --use-std-random --increase-clash-after 1000 --do-recursive"
  cd $OUTDIR
  # Uses OMP
  export OMP_NUM_THREADS=$NPROC
  $SAMPLE_RAMACHANDRAN -i ${samp_dir}/${LIG_SHORT}.cleaned.pdb -o $samp_dir/${LIG_SHORT}_samp -S $hinges_prefix${LIG_SHORT}_fluct_resi.txt $rama_args --do-recursive-fn=$hinges_prefix${LIG_SHORT}_fcc.txt
  $SAMPLE_RAMACHANDRAN -i ${samp_dir}/${REC_SHORT}.cleaned.pdb -o $samp_dir/${REC_SHORT}_samp -S $hinges_prefix${REC_SHORT}_fluct_resi.txt $rama_args --do-recursive-fn=$hinges_prefix${REC_SHORT}_fcc.txt
}

add_sidechain_atoms() {
  export SCWRL4
  ls $samp_dir/${LIG_SHORT}_samp* $samp_dir/${REC_SHORT}_samp* | xargs -n 1 -P $NPROC sh -c '$SCWRL4 -i $1 -o $1.scfix.pdb' sh

  # Need to save these for future steps of the pipeline.
  (cd $samp_dir;
   ls ${LIG_SHORT}_samp*.scfix.pdb > ligands_premin.txt;
   ls ${REC_SHORT}_samp*.scfix.pdb > recepts_premin.txt)
}

amber_minimize() {
  cd $samp_dir
  $AMBERMIN_PAR ligands_premin.txt 1 500 # Just do 1 proc for now.
  $AMBERMIN_PAR recepts_premin.txt 1 500 # The srun command hates more that that.

  for file in `cat ligands_premin.txt`; do
    echo $file $($AMBEREN $file)
  done > ligands_en.txt
  for file in `cat recepts_premin.txt`; do
    echo $file $($AMBEREN $file)
  done > recepts_en.txt
  # Sometimes Amber minimization fails, so let's strip off anything without energy.
  awk '{if (NF > 1) print}' ligands_en.txt | sort -k 2g -o ligands_en.txt
  awk '{if (NF > 1) print}' recepts_en.txt | sort -k 2g -o recepts_en.txt

  # Now need to sort them by energy, then only use the top k, but randomly sort them.
  head ligands_en.txt -n $NUM_SAMPS | sort -R > ligands_use.txt
  head recepts_en.txt -n $NUM_SAMPS | sort -R > recepts_use.txt
}

cleanup_pdbs() {
  # Only need to align proteins we'll actually use.
  IFS=$'\n' LIG_LIST=($(cat $samp_dir/ligands_use.txt | cut -f 1 -d ' ' | sed 's/.pdb$//'))
  IFS=$'\n' REC_LIST=($(cat $samp_dir/recepts_use.txt | cut -f 1 -d ' ' | sed 's/.pdb$//'))

  # Create these if needs be.
  rm -rf $OUTDIR/aligned
  mkdir $OUTDIR/aligned

  # Copy original here.
  $FIX_PDB $RECEPTOR > $OUTDIR/aligned/orig_receptor.pdb
  $FIX_PDB $LIGAND > $OUTDIR/aligned/orig_ligand.pdb

  for i in `seq 0 $(($NUM_SAMPS - 1))`; do
    lig=${LIG_LIST[$i]}_ambermin.pdb
    rec=${REC_LIST[$i]}_ambermin.pdb
    # Not needed
    ##$FIX_PDB $samp_dir/AMBER/$rec > $OUTDIR/aligned/test_rec.pdb
    ##$FIX_PDB $samp_dir/AMBER/$lig > $OUTDIR/aligned/$lig # Don't need to change this one.
    ##$ALIGN_CONTACT $OUTDIR/orig/rec_gold.pdb $OUTDIR/orig/lig_gold.pdb $OUTDIR/aligned/test_rec.pdb $OUTDIR/aligned/$rec
    $FIX_PDB $samp_dir/AMBER/$rec > $OUTDIR/aligned/$rec
    $FIX_PDB $samp_dir/AMBER/$lig > $OUTDIR/aligned/$lig
  done
}

###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Make sure there's not something funky about the BIOMT things...


# First, generate the samples.
one_time=$(date +%s)
echo "#####################################################"
echo "# CG: 1. Generating samples"
echo "#####################################################"
# Some logic doing with the sampling protocol:
if [[ "$PROTOCOL" == "atom" ]]; then
  gen_atom_samples
elif [[ "$PROTOCOL" == "hinge" ]]; then
  gen_hinge_samples
else
  gen_f3d_samples
fi

# Add sidechains
two_time=$(date +%s)
echo "############################################"
echo "# CG: Time for step 1 is `date -u -d @$(($two_time - $one_time)) +"%T"`# "
echo "############################################"
echo
echo
echo "#####################################################"
echo "# CG: 2. Adding sidechains"
echo "#####################################################"
# Only do this if using f3d, since the others keep all the sidechain atoms on.
if [[ "$PROTOCOL" == "f3d" ]]; then
  add_sidechain_atoms
fi

# Minimize, get everything to ligands_en.txt
three_time=$(date +%s)
echo "############################################"
echo "# CG: Time for step 2 is `date -u -d @$(($three_time - $two_time)) +"%T"`# "
echo "############################################"
echo
echo
echo "#####################################################"
echo "# CG: 3. Minimizing with Amber"
echo "#####################################################"
amber_minimize

four_time=$(date +%s)
echo "############################################"
echo "# CG: Time for step 3 is `date -u -d @$(($four_time - $three_time)) +"%T"`# "
echo "############################################"
echo
echo
echo "#####################################################"
echo "# CG: 3. Doing some cleanup."
echo "#####################################################"
cleanup_pdbs
