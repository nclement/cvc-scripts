PDB=$1
SIZE=$2

    AREA=`cat OUT/${PDB}.area`
   
    VOLUME=`cat OUT/${PDB}.volume`
    
    LJ=`grep fastP OUT/${PDB}-ilj.out | sed -e s/fastP\ \=//`

    CP=`grep fastE OUT/${PDB}-icp.out | sed -e s/fastE\ \=//`

    GB=`grep G\_pol OUT/${PDB}-gpol.out | sed -e s/G\_pol\ \=// -e s/kcal.mol//`

    #echo "[${AREA0}]*0.003 + [${VOLUME0}]*0.035 + [${LJ0}] + [${CP0}] + [${GB0}]*100"
    # TOTAL=`echo "${AREA0}*0.003 + \
    #       ${VOLUME0}*0.035 + \
    #       ${LJ0} + \
    #       ${CP0} + \
    #       ${GB0}*100" | bc`
    TOTAL=`echo "${AREA}*0.22 + \
             ${VOLUME}*0.20 + \
             ${LJ}*0.027 + \
             ${CP}*0.094 +  \
             ${GB}*-0.0037" | bc`

    echo ${PDB} ${TOTAL}
