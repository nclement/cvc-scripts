PDBLIST=$1

# Set to 1 for verbose output
VERBOSE=0
while read SINGLE
do
    PDB=(`echo ${SINGLE} | tr " " "\n"`)
    
    AREA0=`cat OUT/${PDB}.area`
   
    VOLUME0=`cat OUT/${PDB}.volume`
    
    LJ0=`grep fastP OUT/${PDB}-ilj.out | sed -e s/fastP\ \=//`

    CP0=`grep fastE OUT/${PDB}-icp.out | sed -e s/fastE\ \=//`

    GB0=`grep G\_pol OUT/${PDB}-gpol.out | sed -e s/G\_pol\ \=// -e s/kcal.mol//`

    if [[ "$VERBOSE" -eq 1 ]]; then
      # echo "[${AREA0}]*0.003 + \
      #       [${VOLUME0}]*0.035 + \
      #       [${LJ0}] + \
      #       [${CP0}] + \
      #       [${GB0}]*100" \| bc
      echo "[${AREA}]*0.22 + \
            [${VOLUME}]*0.20 + \
            [${LJ}]*0.027 + \
            [${CP}]*0.094 +  \
            [${GB}]*-0.0037"
    fi

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

done <${PDBLIST}
