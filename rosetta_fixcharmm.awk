# This file is rosetta_fixcharmm.awk and runs under linux.
# Useage: awk -f rosetta_fixcharmm.awk fixed_charmmout.pdb.
# Charmm29 (Later versions may be better.) "mis" labels H's, e.g. _HD1 rather
# than 1HD_ and causes Rosetta to misread proline as hydroxyproline.
# Relabels H's from Charmm of PRO from, e.g. _HD1 to 1HD_
# Fixes label on atom CD1 of ILE.
# Relabels HSD and HSE from Charmm as HIS
BEGIN {FIELDWIDTHS=" 12 4 1 3 46"}
{
  if ($2 == " HD1" && $4 == "PRO")
    $2 = "1HD "
  if ($2 == " HD2" && $4 == "PRO")
    $2 = "2HD "
  if ($2 == " HB1" && $4 == "PRO")
    $2 = "1HB "
  if ($2 == " HB2" && $4 == "PRO")
    $2 = "2HB "
  if ($2 == " HG1" && $4 == "PRO")
    $2 = "1HG "
  if ($2 == " HG2" && $4 == "PRO")
    $2 = "2HG "
  if ($2 == " CD " && $4 == "ILE")
    $2 = " CD1"
  if (($4 == "HSD") || ($4 == "HSE"))
    $4 = "HIS"
  printf "%-12s", $1
  printf "%-4s", $2
  printf "%-1s", $3
  printf "%-3s", $4
  printf "%-46s\n", $5

}
