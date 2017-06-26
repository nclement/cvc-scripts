from os import path
from sys import argv
my_argv = argv[1:]
# Print usage and quit.
if len(my_argv) < 2:
    sys.exit("usage: pymol createBonds_pymol.py -- <PDB> <HBONDS_FILE>")

pdb = my_argv[0]
hbonds = my_argv[1]

fn, fext = path.splitext(pdb)

#cmd.load(pdb)
print "load " + pdb + ", " + fn
HFILE = open(hbonds, "r")
# Ignore the first line
HFILE.readline()
for line in HFILE:
    if line[0] == "#":
        continue
    els = line.split()
    # from, to, energy
    f = "id " + els[0]
    t = "id " + els[3]
    energy = str(-1*float(els[8]))

    #print "cgo_arrow "+f+", "+ t +", gap=0, radius="+energy+",hradius="+energy
    print "cgo_arrow "+f+", "+ t +", gap=0, radius="+energy+", hradius=0, hlength=0"

print "zoom " + fn
