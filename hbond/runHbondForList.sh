HBOND_SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
SCRIPTS_DIR="$HBOND_SCRIPTS_DIR/../"
source $SCRIPTS_DIR/Makefile.def

PDBLIST=$1

while read PDB
do
	echo "Processing ${PDB}..."
  lc_PDB=$(echo $PDB | perl -ne 'print lc')
  if [ $lc_PDB != $PDB ]; then 
    echo "WARNING: PDB file should be lower case!!!"
    cp $PDB $lc_PDB
    PDB=$lc_PDB
  fi

  # Call this script will all chains.
  $HBOND_SCRIPTS_DIR/runHbondSingle.sh $PDB ?
done < ${PDBLIST}
