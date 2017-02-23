#! /bin/bash

# Will output the sorted (best, i.e. highest score first) list of scores, including
# the "score scale down factor".

F2OUT=$1

ssd=`grep "score scale down factor" $F2OUT | sed 's/.* =//'`
grep "^#" -v $F2OUT | tr -s ' ' | cut -f 3 -d ' ' | awk "{print \$1 * $ssd}" | tac
