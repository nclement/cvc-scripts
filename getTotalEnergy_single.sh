PDB=$1
SIZE=$2

    AREA0=`cat OUT/${PDB}.area`
   
    VOLUME0=`cat OUT/${PDB}.volume`
    
    LJ0=`grep fastP OUT/${PDB}-ilj.out | sed -e s/fastP\ \=//`

    CP0=`grep fastE OUT/${PDB}-icp.out | sed -e s/fastE\ \=//`

    GB0=`grep G\_pol OUT/${PDB}-gpol.out | sed -e s/G\_pol\ \=// -e s/kcal.mol//`

    #echo "[${AREA0}]*0.003 + [${VOLUME0}]*0.035 + [${LJ0}] + [${CP0}] + [${GB0}]*100"
    TOTAL=`echo "${AREA0}*0.003 + \
          ${VOLUME0}*0.035 + \
          ${LJ0} + \
          ${CP0} + \
          ${GB0}*100" | bc`

    echo ${PDB} ${TOTAL}
