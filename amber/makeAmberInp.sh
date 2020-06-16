CYC=${1:-2000}
NCYC=${2:-500}

if [[ $NCYC -gt $CYC ]]; then
  NCYC=$(($CYC/2))
fi

# No idea why this no longer works.
## echo "refine molecule prior to MD"
## echo " &cntrl"
## echo "  imin = 1,"      # minimization run
## echo "  ntx = 1,"       # read coordinates but not velocities from Inpcrd coord file
## echo "  maxcyc = $CYC," # Maximum minimization cycles
## echo "  ncyc = $NCYC,"    # Uses steepest descent algorithm for first ncyc cycles, then
##                         # conj gradient for ncyc-maxcyc cycles
## echo "  ntb = 0,"       # System is not periodic
## echo "  igb = 0,"       # No GB solvation model
## echo "  cut = 12"       # nonbonded cutoff distance in Angstroms (do NOT reduce below 8.0)
## echo "/"
echo "Brief minimization"
echo "&cntrl"
echo "  imin=1,"
echo "  ntx=1,"
echo "  irest=0,"
echo "  maxcyc=$CYC,"
echo "  ncyc=$NCYC,"
echo "  ntpr=100,"
echo "  ntwx=0,"
echo "  cut=8.0,"
echo "/"
