# Needs to be run with pymol
# Program will compute the contact residues, and then print out a file usable
# by F2Dock for _rmsd_atoms.txt, with only the CA atoms.
#
# Parameters are: <lig> <rec> [distance]
#
#
# Example run is:
#   $PYMOL -qrc get_contact_ca_rmsd_atoms.py 1A2K_l_b.pdb 1A2K_r_b.pdb | grep -v "^PyMOL" > 1A2K_rmsd_atoms.txt
#
# 
from pymol import cmd 
from sys import argv

my_argv = argv[1:]
cmd.load(my_argv[0], "lig")
cmd.load(my_argv[1], "rec")
distance = "5"
if len(my_argv) > 2:
  distance = my_argv[2]

cmd.select("contact", "(lig & (all within " + distance + " of rec)) & name ca")
my_dict = { 'pA' : [] }
cmd.iterate_state(1, "contact","pA.append([ID,x,y,z])",space=my_dict)

print(cmd.count_atoms("lig"))
print(len(my_dict['pA']))
for v in my_dict['pA']:
  print("%d %.3f %.3f %.3f" % (v[0], v[1], v[2], v[3]))
