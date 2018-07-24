""" Compute Contact RMSD the hard way.

Pymol and many other programs don't give us an accurate RMSD, even when the
sequences are (nearly) identical. So we're going to do this all manually here.

The following steps will be used:

  1. Compute all the contact RMSD atoms, which is any residue (Ca atom) on one
  chain (here "receptor" and "ligand") that has any atom within <DIST> of any
  other atom on any other residue.

  2. Compute the sequence alignment of the Receptor and Ligand (chains
  separately).

  3. For each contact residue (determined in Step 1 above), compute the
  distance to the aligned residue (determined in Step 2) and report this as the
  contact RMSD.
"""

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SeqUtils import seq1
from contextlib import contextmanager
from multiprocessing import Pool
from math import log
import sys, os


## if len(sys.argv) < 4:
##   print("usage getcRMSD_pymol_full.py <GOLD_R> <TEST_R> <GOLD_L> <TEST_L> [LIMIT=10]")
##   os._exit(1)
## protR=sys.argv[1]  # Receptor
## protRp=sys.argv[2] # sampled receptor
## protL=sys.argv[3]  # Ligand
## protLp=sys.argv[4] # sampled ligand
## # Vreven et al 2015 (zlab5) suggests to use 10\AA as the limit
## X=10
## if len(sys.argv) > 5:
##   X=int(sys.argv[5])

_VERBOSE=0

def getAlignmentList(str):
  alns = []
  match_idx = 0
  for c in str:
    if c != '-':
      alns.append(match_idx)
      match_idx += 1
    else:
      alns.append(-1)

def testContactRes(aln_list, contact_idx, full_seq):
  test_cont = []
  for x in contact_idx:
    test_cont.append(full_seq[x])
  return test_cont

def getSequence(selection):
  myspace = {'residues': []}
  cmd.iterate("%s & n. ca" % selection, "residues.append(resn)", space=myspace)
  return seq1(''.join(myspace['residues']))

def gap_function(x, y):  # x is gap position in seq, y is gap length
  if y == 0:  # No gap
    return 0
  elif y == 1:  # Gap open penalty
      return -2
  return - (2 + y/4.0 + log(y)/2.0)
blosum_matrix = matlist.blosum62
def getBestAlignment(seq1, seq2):
  if _VERBOSE:
    print "Doing the alignment..."
  alns = pairwise2.align.globaldc(seq1, seq2, blosum_matrix,
      gap_function, gap_function, one_alignment_only=True)
  if _VERBOSE:
    print "Finished!"
  return alns[0]

def getContactRes(sel1, sel2, X):
  cmd.select('contact_rec', '(%s) & (all within %s of (%s)) & n. ca'
      %(sel1, X, sel2))
  sp = {'contact': []}
  cmd.iterate('contact_rec', "contact.append((chain, resn, resi, index))", space=sp)
  return sp['contact']

def getFullRes(sel):
  sp = {'full': []}
  cmd.iterate("%s & n. ca" % sel, "full.append((chain, resn, resi, index))",
      space=sp)
  return sp['full']

def getContactResIndx(contact, full):
  contact_idx = [-1]*len(contact)
  for i, e in enumerate(contact):
    #print "Searching for %d" % i, e
    for j, x in enumerate(full):
      if e == x:
        #print "Found %d at %d!" % (i, j)
        contact_idx[i] = j
        break
      #else:
        #print "Skipping %d" % j, x

  return contact_idx

def SimplifySelection(prot, cont):
  str = "%s & n. ca & " % prot
  chain = None
  prev_resi = -2
  printed_prev = False
  for x in cont:
    this_chain = x[0]
    this_resi = int(x[2])
    if not chain:
      str += "(chain %s & resi " % this_chain
    elif chain != this_chain:
      str += ") & (chain %s & resi " % this_chain

    if this_resi > prev_resi + 1:
      if prev_resi == -2:
        str += "%s" % this_resi
      elif not printed_resi:
        str += "-%s+%s" % (prev_resi, this_resi)
      else:
        str += "+%s" % this_resi
      printed_resi = True
    else:
      printed_resi = False

    chain = this_chain
    prev_resi = this_resi

  if not printed_resi:
    str += "-%s" % prev_resi

  return str + ")"

class PyMolAligner:
  def __init__(self, protR, protRp, protL, protLp, X):
    # Load the proteins.
    cmd.load(protR, 'pR')
    cmd.load(protRp, 'pRp')
    cmd.load(protL, 'pL')
    cmd.load(protLp, 'pLp')

    # Need to do some cleanup of these.
    cmd.alter('all', 'segi=""')
    cmd.remove('not (alt ""+A)')
    cmd.alter('all', 'alt=""')

    self._gold_chains_r = cmd.get_chains('pR')
    self._gold_chains_l = cmd.get_chains('pL')
    self._test_chains_r = cmd.get_chains('pRp')
    self._test_chains_l = cmd.get_chains('pLp')

    # Get the sequences.
    self._seqR_gold = getSequence('pR')
    self._seqL_gold = getSequence('pL')
    self._seqR_test = getSequence('pRp')
    self._seqL_test = getSequence('pLp')

    if _VERBOSE:
      print "R_gold is", self._seqR_gold
      print "R_test is", self._seqR_test
    self._R_aln = getBestAlignment(self._seqR_gold, self._seqR_test)
    self._L_aln = getBestAlignment(self._seqL_gold, self._seqL_test)
    if _VERBOSE:
      print "R_aln is", self._R_aln
      print "L_aln is", self._L_aln

    # Need to get contact residues for gold, then corresponding res from test.
    self._cont_R = getContactRes('pR', 'pL', X)
    self._cont_L = getContactRes('pL', 'pR', X)

    self._full_R = getFullRes('pR')
    self._full_L = getFullRes('pL')
    self._full_Rp = getFullRes('pRp')
    self._full_Lp = getFullRes('pLp')

    self._cont_R_idx = getContactResIndx(self._cont_R, self._full_R)
    self._cont_L_idx = getContactResIndx(self._cont_L, self._full_L)

    aln_list_R = getAlignmentList(self._R_aln[0])
    aln_list_L = getAlignmentList(self._L_aln[0])

    self._cont_Rp = testContactRes(aln_list_R, self._cont_R_idx, self._full_Rp)
    self._cont_Lp = testContactRes(aln_list_L, self._cont_L_idx, self._full_Lp)
    if _VERBOSE:
      print "R cont is ", self._cont_R
      print "L cont is ", self._cont_L
      print "Rp cont is ", self._cont_Rp
      print "Lp cont is ", self._cont_Lp

  def RMSDSep(self, lig, rec):
    simp_R = SimplifySelection('pR', self._cont_R)
    simp_Rp = SimplifySelection(rec, self._cont_Rp)
    simp_L = SimplifySelection('pL', self._cont_L)
    simp_Lp = SimplifySelection(lig, self._cont_Lp)
    assert cmd.count_atoms(simp_R) == len(self._cont_R)
    assert cmd.count_atoms(simp_Rp) == len(self._cont_Rp)
    assert cmd.count_atoms(simp_L) == len(self._cont_L)
    assert cmd.count_atoms(simp_Lp) == len(self._cont_Lp)
    return cmd.pair_fit(simp_R, simp_Rp, simp_L, simp_Lp)

##     # This stuff isn't needed anymore.
##     aln_strs = []
##     for i in range(len(self._cont_R)):
##       x = self._cont_R[i]
##       y = self._cont_Rp[i]
##       aln_strs.append(' /%s//%s/%s/CA ' %('pR', x[0], x[2]))
##       aln_strs.append(' /%s//%s/%s/CA ' %(rec, y[0], y[2]))
##     for i in range(len(self._cont_L)):
##       x = self._cont_L[i]
##       y = self._cont_Lp[i]
##       aln_strs.append(' /%s//%s/%s/CA ' %('pL', x[0], x[2]))
##       aln_strs.append(' /%s//%s/%s/CA ' %(lig, y[0], y[2]))
## 
##     #for x in aln_strs:
##     #  print x, cmd.count_atoms(x.replace(',','+'))
##     if _VERBOSE:
##       print "num alns is", len(aln_strs)
##       for x in sorted(aln_strs):
##         print " ", x, cmd.count_atoms(x.replace(',','+'))
##       print "num atoms is", cmd.count_atoms("+".join(aln_strs))
##       print "+".join(aln_strs)
## 
##     return cmd.pair_fit(*aln_strs)

def get_parser():
  """Set up arguments.
  
  Returns argparse object, e.g.:
    args = get_parser().parse_args()
    print args.dist
  """
  from argparse import ArgumentParser
  
  parser = ArgumentParser(description="Contact RMSD")
  parser.add_argument("-R", "--gold_r", dest="gold_r", required=True,
                      help="Gold Receptor chain (or molecule)")
  parser.add_argument("-L", "--gold_l", dest="gold_l", required=True,
                      help="Gold Ligand chain (or molecule)")

  parser.add_argument("-r", "--test_r", dest="test_r", required=True,
                      help="Test Receptor chain (or molecule)")
  parser.add_argument("-l", "--test_l", dest="test_l", required=True,
                      help="Test Ligand chain (or molecule)")

  parser.add_argument("-X", "--dist", dest="dist", default=10, type=int,
                      help="Max distance between contact atoms")

  # Regex's for lig+receptor conformations.
  parser.add_argument("--lig_confs", dest="lig_confs", default=None,
      help="Regex of Ligand confs", required=True)
  parser.add_argument("--rec_confs", dest="rec_confs", default=None,
      help="Regex of Receptor confs if different from --test_r")

  parser.add_argument("-P", dest="nproc", default=1, type=int,
      help="Number of processor cores to use")

  return parser

def getRegexFiles(regex):
  import glob
  return glob.glob(regex)

def CleanName(prot_name):
  return prot_name.replace(".", "_")

def runSingle(x):
  lig = x[1]
  lig_name = CleanName(lig[:-4]) # Remove the ".pdb"
  cmd.load(lig, lig_name)
  rec = x[2]
  rec_name = CleanName(rec[:-4])
  if rec_name != 'pRp':
    cmd.load(rec, rec_name)

  rms = x[0].RMSDSep(lig_name, rec_name)
  #print "Finished, rms is %f" % rms
  return (lig_name, rms)

def main():
  args = get_parser().parse_args()

  alnr = PyMolAligner(args.gold_r, args.test_r, args.gold_l, args.test_l,
      args.dist)

  lig_files = getRegexFiles(args.lig_confs)
  if args.rec_confs:
    rec_files = getRegexFiles(args.rec_confs)
  else:
    rec_files = ['pRp.pdb']*len(lig_files)

  assert len(lig_files) == len(rec_files)

  if not _VERBOSE:
    # Get rid of unwanted output.
    cmd.feedback("disable", "executive", "results")
  p = Pool(args.nproc)
  all_rms = p.map(runSingle, [(alnr, lig_files[i], rec_files[i]) for i in
    range(len(lig_files))])

  for i in range(len(lig_files)):
    print "%s %f" % all_rms[i]

if __name__ == "pymol":
  main()
