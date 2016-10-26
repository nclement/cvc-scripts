PDBLIST=$1

SCRIPTS_DIR="$(dirname -- "$(readlink -f -- "$0")")"
COMPUTE_SCRIPT=${SCRIPTS_DIR}/computeAllEnergy.sh

SIZE=$2
while read PDB
do
    $COMPUTE_SCRIPT ${PDB} ${SIZE}
done <${PDBLIST}
