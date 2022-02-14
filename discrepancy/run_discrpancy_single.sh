#! /bin/bash -x

#SBATCH -n 1                 # total number of mpi tasks requested
#SBATCH -p normal            # queue (partition) -- normal, development, etc.
#SBATCH -t 12:00:00          # run time (hh:mm:ss)
#SBATCH -A P-MRI2

ml cxx11/4.9.1

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")/.."
source ${SCRIPTS_DIR}/Makefile.def

GEN_PRGN=$F3DOCK_DIR/testPRGN
GEN_RANDOM=${TEXMOL_BIN}/testSampleGenerator
GEN_HALTON=${TEXMOL_BIN}/fsu.edu.halton
FSU_TO_MINE=${TEXMOL_BIN}/../../samples/fsu.edu_to_mine.pl
NEXT_P_RAND=${TEXMOL_BIN}/../../samples/next_n_rand
STAR_DISC_Walbar=${TA_ALGS}/TA_improved_bardelta
STAR_DISC_Wal=${TA_ALGS}/TA_improved_delta

# Variables need to be set:
# Needs tpy, one of:
#     tpy=std (using standard c++ rand)
#     tpy=random (using mt19937 rand)
#     tpy=prgn.mt or prgn.nisan for mt19937 or nisan PRG
#     tpy=halton
#     tpy=naive
# Also needs N
# Also needs d
# Also needs ptfile
# Optional:
#   OVERWRITE=false (default true)


wal_trials=5

NPROC=16
#M=1024
M=1001
M_MINUS=$(($M-1))
STAR_M=15000
SAME_RUNS=100
OVERWRITE=${OVERWRITE:-true}

runWal() {
  t_pts=$1;
  # Run these both at the same time!
  $STAR_DISC_Walbar -iter 100000 -trials $wal_trials $t_pts &
  pid=$!
  $STAR_DISC_Wal -iter 100000 -trials $wal_trials $t_pts
  wait $pid
}

pts=$ptfile
out=$ptfile.last

if [ "$tpy" == 'random' ]; then
  if [ ! -f $pts ] || $OVERWRITE ; then
    $GEN_RANDOM $M $d $N 0 $pts &> /dev/null
  fi
elif [ "$tpy" == 'prgn.mt' ]; then
  if [ ! -f $pts ] || $OVERWRITE ; then
    $GEN_PRGN $M $d $N 1 > $pts 2> /dev/null
  fi
elif [ "$tpy" == 'prgn.nisan' ]; then
  if [ ! -f $pts ] || $OVERWRITE ; then
    $GEN_PRGN $M $d $N 4 > $pts 2> /dev/null
  fi
elif [ "$tpy" == 'std' ]; then
  if [ ! -f $pts ] || $OVERWRITE ; then 
    $GEN_RANDOM $M $d $N 2 $pts &> /dev/null
  fi
elif [ "$tpy" == 'naive' ]; then
  if [ ! -f $pts ] || $OVERWRITE ; then
    $GEN_RANDOM $M $d $N 5 $pts &> /dev/null
  fi
else  # halton
  if [ ! -f $pts ] || $OVERWRITE; then
    zeros=$(printf '0%.0s ' $(eval echo {1..$d}))
    ones=$(printf '1%.0s ' $(eval echo {1..$d}))
    rands=$(for i in `seq $d`; do echo -n "$RANDOM "; done)
    hin=$ptfile.hin
    hout=$ptfile.hout
    cat <<EOF > $hin
$d
$N
$RANDOM
$rands
$ones
$($NEXT_P_RAND $d)
EOF
    $GEN_HALTON $hout < $hin &> /dev/null
    $FSU_TO_MINE $M $hout > $pts
    sed -i "1s/^/$M $d $N\n/" $pts
  fi
fi


# Gotta convert the pointset first
head="$d $N reals"
ptsWal=$pts.Wal
echo $head > $ptsWal
tail -n $N $pts | perl -anle "for(@F) {printf \"%f \", \$_/$M_MINUS;}; print \"\"" >> $ptsWal
# Run the program.
runWal $ptsWal
