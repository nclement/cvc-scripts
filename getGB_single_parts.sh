#! /bin/bash

PDB=$1
folder=`dirname $PDB`
pdb=`basename $PDB`
AREA=`cat $folder/OUT/$pdb.area`
VOLUME=`cat $folder/OUT/$pdb.volume`
LJ=`grep fastP $folder/OUT/$pdb-ilj.out | sed 's/fastP =//'`
CP=`grep fastE $folder/OUT/$pdb-icp.out | sed 's/fastE =//'`
DISP=`grep 'G_disp' $folder/OUT/$pdb-disp.out | sed -e 's/.*=//' -e 's/kcal.mol//'`
GB=`grep 'G_pol' $folder/OUT/$pdb-gpol.out | sed -e 's/G_pol =//' -e 's/kcal.mol//'`

# TOTAL=`echo "${AREA}*0.003 + \
# 					   ${VOLUME}*0.035 + \
# 						 ${LJ} + \
# 						 ${CP} +  \
# 						 ${GB}*100" | bc`
# echo "${AREA}*0.003 + ${VOLUME}*0.035 + ${LJ} + ${CP} + ${GB}*100 = $TOTAL"
TOTAL=`echo "${AREA}*0.22 + \
					   ${VOLUME}*0.20 + \
						 ${LJ}*0.027 + \
						 ${CP}*0.094 +  \
						 ${GB}*-0.0037" | bc`
echo "${AREA}*0.22 + ${VOLUME}*0.20 + ${LJ}*0.027 + ${CP}0.094 + ${GB}*-0.037 = $TOTAL"
