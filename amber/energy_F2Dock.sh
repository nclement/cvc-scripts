#! /bin/bash

# Computes the energy on the top k predictions from F2dock.

# Just needs F2Dock file and top number of predictions.
F2Dock=$1
TOP_K=$2

# Some programs to use.
EXTRACT_X=/work/01872/nclement/scripts/docking/postprocessing/extractXForms.pl
CONF_GEN=/work/01872/nclement/scripts/docking/postprocessing/conformationGenerator
AMBER_RUN=/work/01872/nclement/scripts/amber/get_amberen_all.sh

# Get the static molecule
static=`grep staticMoleculeF2d $F2Dock | sed 's/.* = \(.*\)_1.7.f2d/\1/'`
moving=`grep movingMoleculeF2d $F2Dock | sed 's/.* = \(.*\).f2d/\1/'`
$EXTRACT_X $TOP_K $F2Dock > xforms.txt
$CONF_GEN $moving.pdb xforms.txt 0 $TOP_K

# Then compute the energy on all of them.
$AMBER_RUN
