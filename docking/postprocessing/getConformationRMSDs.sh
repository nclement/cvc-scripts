#! /bin/bash

POST_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR=$POST_SCRIPTS_DIR/../../
source $SCRIPTS_DIR/Makefile.def

CRMSD=$SCRIPTS_DIR/getcRMSD_pymol_norot.sh

################################################################################
# Often, we'd like to know the energy of the top k docked conformations. This  #
# takes a few steps to complete. Do them all here at once.                     #
################################################################################

REC_GOLD=$1   # Receptor Gold
LIG_GOLD=$2   # Ligand Gold
F2OUT=$3 # F2Dock output file
NUM=$4   # Number of conformations to generate.
OUTFILE=$5 # Output file for final energy.
NPROC=${6:-16} # Number of cores to use.

# Get the running rec+lig from the F2Dock outfile.
REC=$(grep "_1.7.f2d" $F2OUT | sed 's/.* = \(.*\)_1.7.f2d/\1.pdb/' ); \
LIG=$(grep "moving.*f2d" $F2OUT | sed 's/.* = \(.*\).f2d/\1.pdb/' ); \

echo "Running with:"
echo "  REC: $REC GOLD: $REC_GOLD"
echo "  LIG: $LIG GOLD: $LIG_GOLD"
echo "  F2OUT: $F2OUT"
echo "  NUM: $NUM"
echo "  OUTFILE: $OUTFILE"
echo "  NPROC: $NPROC"

WORKDIR=$(readlink -f -- "$OUTFILE")_$(basename $LIG)_$(basename $REC)
echo "Now using [$WORKDIR] as the working directory."
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
#for file in ${LIG_PRE}_conf_*.pdb; do cat $REC >> $file; done;

ls ${LIG::${#LIG}-4}_conf_*.pdb > $tmpSingles;

export CRMSD REC_GOLD REC LIG_GOLD
xargs -a $tmpSingles -P$NPROC -Ipdb sh -c 'echo pdb $($CRMSD $REC_GOLD $REC $LIG_GOLD pdb)' | tee $OUTFILE

# Clean up
#rm -rf $WORKDIR
