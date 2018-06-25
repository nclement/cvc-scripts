#! /bin/bash

POST_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR=$POST_SCRIPTS_DIR/../../
source $SCRIPTS_DIR/Makefile.def


################################################################################
# Often, we'd like to know the energy of the top k docked conformations. This  #
# takes a few steps to complete. Do them all here at once.                     #
################################################################################

REC=$1   # Receptor
LIG=$2   # Ligand
F2OUT=$3 # F2Dock output file
NUM=$4   # Number of conformations to generate.
OUTFILE=$5 # Output file for final energy.
NPROC=${6:-16} # Number of cores to use.
SIZE=256

echo "Running with:"
echo "  REC: $REC"
echo "  LIG: $LIG"
echo "  F2OUT: $F2OUT"
echo "  NUM: $NUM"
echo "  OUTFILE: $OUTFILE"
echo "  NPROC: $NPROC"
echo "  SIZE: $SIZE"

WORKDIR=$(readlink -f -- "$OUTFILE")_$(basename $LIG)_$(basename $REC)
echo "Now using [$WORKDIR] as the working directory."
rm -rf $WORKDIR
mkdir -p $WORKDIR

# Need to get the full path on output
OUTFILE=$(readlink -f -- "$OUTFILE")
# Also need to have these files in the tmp dir.
cp $REC $LIG $WORKDIR
REC=$(basename $REC)
LIG=$(basename $LIG)
# Some tmp files we need.
tmpXforms=$(mktemp $WORKDIR/xforms-tmp.XXXXXX)
tmpSingles=$(mktemp $WORKDIR/singles.XXXXXX)

# First, get the xforms file.
$POST_SCRIPTS_DIR/extractXForms.pl $NUM $F2OUT > $tmpXforms
# Then, generate the conformations (do it in the $WORKDIR dir so we don't pollute the workspace).
cd $WORKDIR;

# They'll all be called ${LIG::${#LIG}-4}_conf_XX.pdb
LIG_PRE=${LIG::${#LIG}-4}
rm ${LIG_PRE}_conf_*
$POST_SCRIPTS_DIR/conformationGenerator $LIG $tmpXforms 0 $(($NUM-1));
# Add the receptor to all of them.
for file in ${LIG_PRE}_conf_*.pdb; do 
  # Try and add a TER at the end of this one if we need to.
  if ! grep -q TER $file; then
    echo "TER" >> $file;
  fi
  cat $REC >> $file;
  # Don't like END in the middle.
  sed -i "s/^END/TER/" $file 
  # Also don't like blank lines
  sed -i '/^$/d' $file
done;

ls ${LIG::${#LIG}-4}_conf_*.pdb > $tmpSingles;

# Get the energy for all the files.
$SCRIPTS_DIR/runSingles_parallel.sh $tmpSingles $SIZE $NPROC
$SCRIPTS_DIR/getTotalEnergy.sh $tmpSingles > $OUTFILE


# Or do this by amber
#echo Running amber on all samples.
#$SCRIPTS_DIR/amber/get_amberen_all.sh "${LIG_PRE}_conf_*.pdb" | grep "^$LIG_PRE" > $OUTFILE.amber

# Clean up
#rm -rf $WORKDIR
