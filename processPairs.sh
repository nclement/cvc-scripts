PDBLIST=$1

SIZE=$2
while read PAIR
do
    PDB=(`echo ${PAIR} | tr " " "\n"`)

    AREA0=`cat OUT/${PDB[0]}.area`
    AREA1=`cat OUT/${PDB[1]}.area`
    AREA01=`cat OUT/${PDB[0]}-${PDB[1]}.area`
   
    VOLUME0=`cat OUT/${PDB[0]}.volume`
    VOLUME1=`cat OUT/${PDB[1]}.volume`
    VOLUME01=`cat OUT/${PDB[0]}-${PDB[1]}.volume`       
    
    LJ0=`grep fastP OUT/${PDB[0]}-ilj.out | sed -e s/fastP\ \=//`
    LJ1=`grep fastP OUT/${PDB[1]}-ilj.out | sed -e s/fastP\ \=//`
    LJ01=`grep fastP OUT/${PDB[0]}-${PDB[1]}-ilj.out | sed -e s/fastP\ \=//`

    #LJ01a=`grep fastP OUT/${PDB[0]}-${PDB[1]}-lj.out | sed -e s/fastP\ \=//`

    #echo ${LJ0} ${LJ1} ${LJ01} ${LJ01a}
    #echo ${LJ01} - ${LJ0} -${LJ1} | bc

    CP0=`grep fastE OUT/${PDB[0]}-icp.out | sed -e s/fastE\ \=//`
    CP1=`grep fastE OUT/${PDB[1]}-icp.out | sed -e s/fastE\ \=//`
    CP01=`grep fastE OUT/${PDB[0]}-${PDB[1]}-icp.out | sed -e s/fastE\ \=//`

    #CP01a=`grep fastE OUT/${PDB[0]}-${PDB[1]}-cp.out | sed -e s/fastE\ \=//`
    #CP01b=`grep naiveE OUT/${PDB[0]}-${PDB[1]}-cpe.out | sed -e s/naiveE\ \=//`
    
    #echo ${CP0} ${CP1} ${CP01} ${CP01a} ${CP01b}
    #echo ${CP01} - ${CP0} -${CP1} | bc

    GB0=`grep G\_pol OUT/${PDB[0]}-gpol.out | sed -e s/G\_pol\ \=// -e s/kcal.mol//`
    GB1=`grep G\_pol OUT/${PDB[1]}-gpol.out | sed -e s/G\_pol\ \=// -e s/kcal.mol//`
    GB01=`grep G\_pol OUT/${PDB[0]}-${PDB[1]}-gpol.out | sed -e s/G\_pol\ \=// -e s/kcal.mol//`
    
    #echo ${AREA0} ${AREA1} ${AREA01}
    #echo ${VOLUME0} ${VOLUME1} ${VOLUME01}
    #echo ${LJ0} ${LJ1} ${LJ01}
    #echo ${CP0} ${CP1} ${CP01}
    #echo ${GB0} ${GB1} ${GB01}

    TOTAL=`echo "(${AREA01} - ${AREA0} - ${AREA1})*0.003 + \
          (${VOLUME01} - ${VOLUME0} - ${VOLUME1})*0.035 + \
         (${LJ01} - ${LJ0} - ${LJ1}) + \
          (${CP01} - ${CP0} - ${CP1}) + \
         (${GB01} - ${GB0} - ${GB1})*100" | bc`

     #echo "(${AREA01} - ${AREA0} - ${AREA1})*0.003 + \
     #     (${VOLUME01} - ${VOLUME0} - ${VOLUME1})*0.035 + \
     #     (${LJ01} - ${LJ0} - ${LJ1}) + \
     #     (${CP01} - ${CP0} - ${CP1}) + \
     #     (${GB01} - ${GB0} - ${GB1})*100"

    #TOTAL=`echo "(${AREA01} - ${AREA0} - ${AREA1})*0.003 + \
    #      (${VOLUME01} - ${VOLUME0} - ${VOLUME1})*0.035 + \
    #      (${LJ01}) + \
    #      (${CP01})*10000 + \
    #      (${GB01} - ${GB0} - ${GB1})*100" | bc`

    #TOTAL=`echo "(${AREA01} - ${AREA0} - ${AREA1})*0.003 + \
    #      (${VOLUME01} - ${VOLUME0} - ${VOLUME1})*0.035 + \
    #      (${GB01} - ${GB0} - ${GB1})*100" | bc`

    #TOTAL=`echo "(${GB01} - ${GB0} - ${GB1})" | bc`



    echo ${PDB[0]} ${PDB[1]} ${TOTAL}

    #cat PQR/${PDBS[0]}.pqr PQR/${PDBS[1]}.pqr >PQR/${PDBS[0]}-${PDBS[1]}.pqr

    #./SH/computeAllEnergy.sh ${PDBS[0]}-${PDBS[1]} ${SIZE}

done <${PDBLIST}
