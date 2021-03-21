#!/bin/bash

SCRIPTS="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS=$SCRIPTS/..
source $SCRIPTS/Makefile.def

PDB=$1        # Full path to PDB
NSAMP=$2      # Num samples to generate
OUT_PREFIX=$3 # Prefix for samples (local or absolute)

TEMP_FACTOR_BINARY=$MOLSURF_BIN/runTempFactorUQ
ADD_RESIDUES=$SCRIPTS/certifyResidueNum.pl

# Generate samples.
echo $TEMP_FACTOR_BINARY $PDB $NSAMP $OUT_PREFIX
$TEMP_FACTOR_BINARY $PDB $NSAMP $OUT_PREFIX

# Sometimes the samples have errors in the output. This should fix them.
for file in $OUT_PREFIX*.pdb; do
  perl -pi -e 's/(\d)-(\d\d\d\.\d\d)\d/$1 -$2/' $file;
done

# Sometimes the files are missing the residue numbers, and this wrecks havoc
# with Amber. Fix them with this script:
for file in $OUT_PREFIX*.pdb; do
  $ADD_RESIDUES $PDB $file $file; # Overwrite the same file.
done

