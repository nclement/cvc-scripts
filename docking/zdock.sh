#! /bin/bash

DOCKING_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR=${DOCKING_SCRIPTS_DIR}/../
source ${SCRIPTS_DIR}/Makefile.def

REC=$1 # Need to run mark_sur on the receptor already or it won't work
LIG=$2
out=`basename ${REC}`_`basename ${LIG}`_zdock.out
rec=${REC%pdb}zdock.pdb
lig=${LIG%pdb}zdock.pdb

cp ${ZDOCK_DIR}/uniCHARMM .
#${ZDOCK_DIR}/mark_sur $REC $rec
${ZDOCK_DIR}/mark_sur $LIG $lig
${ZDOCK_DIR}/zdock -R $rec -L $lig -o $out
#head -n 15 $out | tail -n 10 | tr -s " " | cut -f 7 > $out.score
