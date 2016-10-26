PDBA=$1
PDBB=$2
colorA=$3
colorB=$4
colorInterfaceA=$5
colorInterfaceB=$6

echo "Running for:"
echo "  $PDBA: $colorA -- $colorInterfaceA"
echo "  $PDBB: $colorB -- $colorInterfaceB"

PDBA_noext=${PDBA%.pdb}
PDBB_noext=${PDBB%.pdb}

mkdir -p STATS
mkdir -p RAWNC
mkdir -p PBOUT

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
#SCRIPTS=/work/01872/nclement/scripts/
clashOctree=${SCRIPTS_DIR}/binaries/clashOctree
rescontOctree=${SCRIPTS_DIR}/binaries/rescontOctree
interfaceStats=${SCRIPTS_DIR}/binaries/testInterfaceStats
dummyXforms=$SCRIPTS_DIR/stats/xforms.txt
color=${SCRIPTS_DIR}/binaries/colorify
#MolSurf=/work/01872/nclement/MolSurf/release/bin/MolSurf
getAll=$SCRIPTS_DIR/getTotalEnergy_single.sh

#Read all the energy terms which were computed by computeAllEnergies.
#The D_ variables are the computed change in the energy term.
AREA0=`cat OUT/${PDBA}.area`
AREA1=`cat OUT/${PDBB}.area`
AREA01=`cat OUT/${PDBA}-${PDBB}.area`
DAREA=`echo ${AREA01} - ${AREA1} - ${AREA0} | bc`   

VOLUME0=`cat OUT/${PDBA}.volume`
VOLUME1=`cat OUT/${PDBB}.volume`
VOLUME01=`cat OUT/${PDBA}-${PDBB}.volume`       
DVOLUME=`echo ${VOLUME01} - ${VOLUME1} - ${VOLUME0} | bc`

LJ0=`grep fastP OUT/${PDBA}-ilj.out | sed -e s/fastP\ \=//`
LJ1=`grep fastP OUT/${PDBB}-ilj.out | sed -e s/fastP\ \=//`
LJ01=`grep fastP OUT/${PDBA}-${PDBB}-ilj.out | sed -e s/fastP\ \=//`
DLJ=`echo ${LJ01} - ${LJ1} - ${LJ0} | bc`

CP0=`grep fastE OUT/${PDBA}-icp.out | sed -e s/fastE\ \=//`
CP1=`grep fastE OUT/${PDBB}-icp.out | sed -e s/fastE\ \=//`
CP01=`grep fastE OUT/${PDBA}-${PDBB}-icp.out | sed -e s/fastE\ \=//`
DCP=`echo ${CP01} - ${CP1} - ${CP0} | bc`

DP0=`grep 'G_disp' OUT/${PDBA}-disp.out | sed -e 's/.*=//' -e 's/kcal.mol//'`
DP1=`grep 'G_disp' OUT/${PDBB}-disp.out | sed -e 's/.*=//' -e 's/kcal.mol//'`
DP01=`grep 'G_disp' OUT/${PDBA}-${PDBB}-disp.out | sed -e 's/.*=//' -e 's/kcal.mol//'`
DDISP=`echo $DP01 - $DP1 - $DP0 | bc`

GB0=`grep G\_pol OUT/${PDBA}-gpol.out | sed -e s/G\_pol\ \=// -e s/kcal.mol//`
GB1=`grep G\_pol OUT/${PDBB}-gpol.out | sed -e s/G\_pol\ \=// -e s/kcal.mol//`
GB01=`grep G\_pol OUT/${PDBA}-${PDBB}-gpol.out | sed -e s/G\_pol\ \=// -e s/kcal.mol//`
DGB=`echo ${GB01} - ${GB1} - ${GB0} | bc`

pbdir=OUT
if [ -f OUT/${PDBA}.energy ]; then
  pbdir=OUT
elif [ -f PBOUT/${PDBA}.energy ]; then
  pbdir=PBOUT
elif [ -f PBOUT256/${PDBA}.energy ]; then
  pbdir=PBOUT256
fi

PB0=`grep 'Total Energy' ${pbdir}/${PDBA}.energy | sed -e 's/.*\. //'`
PB1=`grep 'Total Energy' ${pbdir}/${PDBB}.energy | sed -e 's/.*\. //'`
PB01=`grep 'Total Energy' ${pbdir}/${PDBA}-${PDBB}.energy | sed -e 's/.*\. //'`
DPB=`echo $PB01 - $PB1 - $PB0 | bc`

TOTAL0=`${getAll} ${PDBA} | sed -e 's/^.* //'`
TOTAL1=`${getAll} ${PDBB} | sed -e 's/^.* //'`
TOTAL01=`${getAll} ${PDBA}-${PDBB} | sed -e 's/^.* //'`
DTOTAL=`echo $TOTAL01 - $TOTAL1 - $TOTAL0 | bc`

#Output the data in a nice tex-friendly format
echo "---Computing energy terms---"
echo "\\begin{table}"
echo "\\begin{tabular}{l|rrrr}"
echo "Energy Term & Combined & ${PDBA} & ${PDBB} & Change \\\\"
echo "AREA & $AREA01 & $AREA0 & $AREA1 & $DAREA \\\\"
echo "VOLUME & $VOLUME01 & $VOLUME0 & $VOLUME1 & $DVOLUME \\\\"
echo "LJ & $LJ01 & $LJ0 & $LJ1 & $DLJ \\\\"
echo "CP & $CP01 & $CP0 & $CP1 & $DCP \\\\"
echo "DISP & $DP01 & $DP0 & $DP1 & $DDISP \\\\"
echo "GB & $GB01 & $GB0 & $GB1 & $DGB \\\\"
echo "PB & $PB01 & $PB0 & $PB1 & $DPB \\\\"
echo "Total G & $TOTAL01 & $TOTAL0 & $TOTAL1 & $DTOTAL \\\\"
echo "\\end{tabular}"
echo "\\end{table}"
echo ""

#Compute clashes, output as output goes
echo "---Clashes---"
echo `${clashOctree} PQR/${PDBA}.pqr PQR/${PDBB}.pqr`
echo ""
#Compute resides, output as output goes
echo "---Residue-Residue contact score---"
echo `${rescontOctree} PQR/${PDBA}.pqr PQR/${PDBB}.pqr`
echo ""

echo "---Interface Stats---"
#Run interface stats and store in a temp file
echo ${interfaceStats} PQR/${PDBA}.pqr PQR/${PDBB}.pqr QUAD/${PDBA}.quad QUAD/${PDBB}.quad ${dummyXforms} \> STATS/tempstats
${interfaceStats} PQR/${PDBA}.pqr PQR/${PDBB}.pqr QUAD/${PDBA}.quad QUAD/${PDBB}.quad ${dummyXforms} > STATS/tempstats
${interfaceStats} PQR/${PDBA}.pqr PQR/${PDBB}.pqr QUAD/${PDBA}.quad QUAD/${PDBB}.quad ${dummyXforms} > STATS/tempstats
#Change colons to &, for tex format
tr ':' '&' < STATS/tempstats > STATS/tempstats2
#Replace newlines with \\ for tex format
awk 1 ORS='\\\\\n' STATS/tempstats2 > STATS/tempstats3
#Keep only the lines we're interested in 
sed -n '3,5p' STATS/tempstats3 > STATS/tempstats4
sed -n '14,19p' STATS/tempstats3 >> STATS/tempstats4
#Display
cat STATS/tempstats4
echo ""

#Here we color the surfaces, and compute the interface surface.
echo "---Coloring Surfaces and computing interfaces. First surface color ${colorA}.  Second surface color ${colorB}.---"
echo "---Output in RAWNC/---"
echo $color RAWN/${PDBA}.rawn RAWNC/${PDBA}.rawnc $colorA 
$color RAWN/${PDBA}.rawn RAWNC/${PDBA}.rawnc $colorA 
echo $color RAWN/${PDBB}.rawn RAWNC/${PDBB}.rawnc $colorB 
$color RAWN/${PDBB}.rawn RAWNC/${PDBB}.rawnc $colorB 
echo $MolSurf -getInterfaceSurface RAWNC/${PDBA}.rawnc RAWNC/${PDBB}.rawnc 4.0 1 $colorInterfaceA $colorInterfaceB
$MolSurf -getInterfaceSurface RAWNC/${PDBA}.rawnc RAWNC/${PDBB}.rawnc 4.0 1 $colorInterfaceA $colorInterfaceB
mv RAWNC/${PDBA_noext}_interface.rawnc RAWNC/${PDBA_noext}-${PDBB_noext}_interface.rawnc
mv RAWNC/${PDBB_noext}_interface.rawnc RAWNC/${PDBB_noext}-${PDBA_noext}_interface.rawnc
$MolSurf -getInterfaceSurface PBOUT/$PDBA.rawnc PBOUT/$PDBB.rawnc 4.0 1
mv PBOUT/${PDBA_noext}_interface.rawnc PBOUT/${PDBA_noext}-${PDBB_noext}_interface.rawnc
mv PBOUT/${PDBB_noext}_interface.rawnc PBOUT/${PDBB_noext}-${PDBA_noext}_interface.rawnc
echo "---Now you can open up RAWNC/*_*-interface.rawnc and PBOUT/*_*-interface.rawnc
and take some pictures.---"
