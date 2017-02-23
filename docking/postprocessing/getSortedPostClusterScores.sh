#! /bin/bash

# Will output the sorted (best, i.e. highest score first) list of scores, including
# the "score scale down factor".
# This script uses the "score after cluster" from the previous version of F2dock.

F2OUT=$1 # Output for F2Dock
if !  grep -q "COLNAME score after cluster" $F2OUT ; then
  echo "Error: Probably using the wrong F2Dock version!!"
  echo "    cannot find score after cluster"
else
  ssd=`grep "score scale down factor" $F2OUT | sed 's/.* =//'`
  grep "^#" -v $F2OUT | tr -s ' ' | cut -f 9 -d ' ' | awk "{print \$1 * $ssd}" | tac
fi
