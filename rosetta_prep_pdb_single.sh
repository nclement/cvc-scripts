#! /bin/bash

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
ADD_TER=$SCRIPTS_DIR/pdb_add_ter.pl

PDB=$1
out=$2
rl=$3

grep "^\(ATOM\|HETATM\|TER\)" $PDB | sed 's/^\(.\{21\}\)./\1'${rl}'/' | $ADD_TER > $out
awk -f $SCRIPTS_DIR/rosetta_fixcharmm.awk $out > ${out%.pdb}_cleaned.pdb
mv ${out%.pdb}_cleaned.pdb $out
