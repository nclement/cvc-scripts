#! /bin/bash

# Input is the PDB file entered to runAmber_single.sh
# The output to this is located in something like:
#  AMBER/${PDB}.amberout

PDB=$1

AMBEROUT=AMBER/$PDB.amberout

grep -A 20 "FINAL RESULTS" $AMBEROUT | grep -A 1 ENERGY | tail -n 1 | tr -s ' ' | cut -f 3 -d ' '
