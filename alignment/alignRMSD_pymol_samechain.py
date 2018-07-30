""" Compute Contact RMSD the hard way.

Pymol and many other programs don't give us an accurate RMSD, even when the
sequences are (nearly) identical. So we're going to do this all manually here.

The following steps will be used:

  1. Compute all the contact RMSD atoms, which is any residue (Ca atom) on one
  chain (here 'receptor' and 'ligand') that has any atom within <DIST> of any
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
_REC_TOGETHER="__tog"

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
  cmd.iterate('%s & n. ca' % selection, 'residues.append(resn)', space=myspace)
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
    print 'Doing the alignment...'
  alns = pairwise2.align.globaldc(seq1, seq2, blosum_matrix,
      gap_function, gap_function, one_alignment_only=True)
  if _VERBOSE:
    print 'Finished!'
  return alns[0]

def getContactRes(sel1, sel2, X):
  cmd.select('contact_rec', '(%s) & (all within %s of (%s)) & n. ca'
      %(sel1, X, sel2))
  sp = {'contact': []}
  cmd.iterate('contact_rec', 'contact.append((chain, resn, resi, index))', space=sp)
  return sp['contact']

def getFullRes(sel):
  sp = {'full': []}
  cmd.iterate('%s & n. ca' % sel, 'full.append((chain, resn, resi, index))',
      space=sp)
  return sp['full']

def getContactResIndx(contact, full):
  contact_idx = [-1]*len(contact)
  for i, e in enumerate(contact):
    #print 'Searching for %d' % i, e
    for j, x in enumerate(full):
      if e == x:
        #print 'Found %d at %d!' % (i, j)
        contact_idx[i] = j
        break
      #else:
        #print 'Skipping %d' % j, x

  return contact_idx

def SimplifySelection(prot, cont):
  str = '%s & n. ca & ' % prot
  chain = None
  prev_resi = -2
  printed_prev = False
  for x in cont:
    this_chain = x[0]
    this_resi = int(x[2])
    if not chain:
      str += '(chain %s & resi ' % this_chain
    elif chain != this_chain:
      str += ') & (chain %s & resi ' % this_chain

    if this_resi > prev_resi + 1:
      if prev_resi == -2:
        str += '%s' % this_resi
      elif not printed_resi:
        str += '-%s+%s' % (prev_resi, this_resi)
      else:
        str += '+%s' % this_resi
      printed_resi = True
    else:
      printed_resi = False

    chain = this_chain
    prev_resi = this_resi

  if not printed_resi:
    str += '-%s' % prev_resi

  return str + ')'

class PyMolAligner:
  def __init__(self, protR, protRp, protL, protLp, X):
    # Remember these.
    self._protR = protR
    self._protL = protL
    self._protRp = protRp
    self._protLp = protLp

    self._gold_chains_r = cmd.get_chains(protR)
    self._gold_chains_l = cmd.get_chains(protL)
    self._test_chains_r = cmd.get_chains(protRp)
    self._test_chains_l = cmd.get_chains(protLp)

    # Get the sequences.
    self._seqR_gold = getSequence(protR)
    self._seqL_gold = getSequence(protL)
    self._seqR_test = getSequence(protRp)
    self._seqL_test = getSequence(protLp)

    if _VERBOSE:
      print 'R_gold is', self._seqR_gold
      print 'R_test is', self._seqR_test
    self._R_aln = getBestAlignment(self._seqR_gold, self._seqR_test)
    self._L_aln = getBestAlignment(self._seqL_gold, self._seqL_test)
    if _VERBOSE:
      print 'R_aln is', self._R_aln
      print 'L_aln is', self._L_aln

    # Need to get contact residues for gold, then corresponding res from test.
    self._cont_R = getContactRes(protR, protL, X)
    self._cont_L = getContactRes(protL, protR, X)

    self._full_R = getFullRes(protR)
    self._full_L = getFullRes(protL)
    self._full_Rp = getFullRes(protRp)
    self._full_Lp = getFullRes(protLp)

    self._cont_R_idx = getContactResIndx(self._cont_R, self._full_R)
    self._cont_L_idx = getContactResIndx(self._cont_L, self._full_L)

    aln_list_R = getAlignmentList(self._R_aln[0])
    aln_list_L = getAlignmentList(self._L_aln[0])

    self._cont_Rp = testContactRes(aln_list_R, self._cont_R_idx, self._full_Rp)
    self._cont_Lp = testContactRes(aln_list_L, self._cont_L_idx, self._full_Lp)
    if _VERBOSE:
      print 'R cont (len:%d) is ' % len(self._cont_R), self._cont_R
      print 'L cont (len:%d) is ' % len(self._cont_L), self._cont_L
      print 'Rp cont (len:%d) is ' % len(self._cont_Rp), self._cont_Rp
      print 'Lp cont (len:%d) is ' % len(self._cont_Lp), self._cont_Lp

  def RMSDTog(self, prot):
    simp_R = SimplifySelection(self._protR, self._cont_R)
    simp_Rp = SimplifySelection(prot, self._cont_Rp)
    simp_L = SimplifySelection(self._protL, self._cont_L)
    simp_Lp = SimplifySelection(prot, self._cont_Lp)
    if _VERBOSE:
      print '%s: %d vs %d' %(simp_R, cmd.count_atoms(simp_R), len(self._cont_R))
      print '%s: %d vs %d' %(simp_Rp, cmd.count_atoms(simp_Rp), len(self._cont_Rp))
      print '%s: %d vs %d' %(simp_L, cmd.count_atoms(simp_L), len(self._cont_L))
      print '%s: %d vs %d' %(simp_Lp, cmd.count_atoms(simp_Lp), len(self._cont_Lp))
    assert cmd.count_atoms(simp_R) == len(self._cont_R)
    assert cmd.count_atoms(simp_Rp) == len(self._cont_Rp)
    assert cmd.count_atoms(simp_L) == len(self._cont_L)
    assert cmd.count_atoms(simp_Lp) == len(self._cont_Lp)
    return cmd.pair_fit(simp_R, simp_Rp, simp_L, simp_Lp)

  def RMSDSep(self, rec, lig):
    simp_R = SimplifySelection(self._protR, self._cont_R)
    simp_Rp = SimplifySelection(rec, self._cont_Rp)
    simp_L = SimplifySelection(self._protL, self._cont_L)
    simp_Lp = SimplifySelection(lig, self._cont_Lp)
    if _VERBOSE:
      print '%s: %d vs %d' %(simp_R, cmd.count_atoms(simp_R), len(self._cont_R))
      print '%s: %d vs %d' %(simp_Rp, cmd.count_atoms(simp_Rp), len(self._cont_Rp))
      print '%s: %d vs %d' %(simp_L, cmd.count_atoms(simp_L), len(self._cont_L))
      print '%s: %d vs %d' %(simp_Lp, cmd.count_atoms(simp_Lp), len(self._cont_Lp))
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
##       aln_strs.append(' /%s//%s/%s/CA ' %(self._protR, x[0], x[2]))
##       aln_strs.append(' /%s//%s/%s/CA ' %(rec, y[0], y[2]))
##     for i in range(len(self._cont_L)):
##       x = self._cont_L[i]
##       y = self._cont_Lp[i]
##       aln_strs.append(' /%s//%s/%s/CA ' %(self._protL, x[0], x[2]))
##       aln_strs.append(' /%s//%s/%s/CA ' %(lig, y[0], y[2]))
## 
##     #for x in aln_strs:
##     #  print x, cmd.count_atoms(x.replace(',','+'))
##     if _VERBOSE:
##       print 'num alns is', len(aln_strs)
##       for x in sorted(aln_strs):
##         print ' ', x, cmd.count_atoms(x.replace(',','+'))
##       print 'num atoms is', cmd.count_atoms('+'.join(aln_strs))
##       print '+'.join(aln_strs)
## 
##     return cmd.pair_fit(*aln_strs)

def setSplitArgs(parser):
  parser.add_argument('-R', '--gold_r', dest='gold_r', required=True,
                      help='Gold Receptor chain (or molecule)')
  parser.add_argument('-L', '--gold_l', dest='gold_l', required=True,
                      help='Gold Ligand chain (or molecule)')

  parser.add_argument('-r', '--test_r', dest='test_r', required=True,
                      help='Test Receptor chain (or molecule)')
  parser.add_argument('-l', '--test_l', dest='test_l', required=True,
                      help='Test Ligand chain (or molecule)')

  # Regex's for lig+receptor conformations.
  parser.add_argument('--lig_confs', dest='lig_confs', default=None,
      help='Regex of Ligand confs')
  parser.add_argument('--rec_confs', dest='rec_confs', default=None,
      help='Regex of Receptor confs if different from --test_r')
  # Also possible to specify only a single one. But this makes things tricky.
  parser.add_argument('--prot_confs', dest='prot_confs', default=None,
      help='Regex of proteins if specified together '+
           '(ignores rec_confs and lig_confs).')

def setTogArgs(parser):
  # For when the PDBs are not split.
  parser.add_argument('-P', '--gold', dest='gold_p', required=True,
          help='Gold protein to test against')
  parser.add_argument('-p', '--test', dest='test_p', required=True,
          help='Model test protein')

  # Also need chains.
  parser.add_argument('-R', '--gold_rec_chains', dest='gold_rec_chains',
          required=True, help='Chains of gold receptor')
  parser.add_argument('-L', '--gold_lig_chains', dest='gold_lig_chains',
          required=True, help='Chains of gold ligand')
  parser.add_argument('-r', '--test_rec_chains', dest='test_rec_chains',
          required=True, help='Chains of test receptor')
  parser.add_argument('-l', '--test_lig_chains', dest='test_lig_chains',
          required=True, help='Chains of test ligand')

  # Input prots.
  parser.add_argument('--confs', dest='test_confs', default=None, required=True,
          help='Regex of test protein configurations')


def addCommonArgs(parser):
  # Vreven et al 2015 (zlab5) suggests to use 10\AA as the limit
  parser.add_argument('-X', '--dist', dest='dist', default=10, type=int,
                      help='Max distance between contact atoms')
  parser.add_argument('-N', dest='nproc', default=1, type=int,
      help='Number of processor cores to use')
  parser.add_argument('--input_fn', dest='input_fn', default=None,
      help='Use if inputs to --*_confs are actually a filename, not a Regex',
      action='store_true')

def get_parser():
  """Set up arguments.
  
  Returns argparse object, e.g.:
    args = get_parser().parse_args()
    print args.dist
  """
  from argparse import ArgumentParser
  
  parser = ArgumentParser(description='Contact RMSD')

  subparsers = parser.add_subparsers(help='PDB input options')
  parser_split = subparsers.add_parser(
          'split_pdb', help='Use when input PDBs are split')
  setSplitArgs(parser_split)
  parser_tog = subparsers.add_parser(
          'tog_pdb', help='Use when input PDBs are in the same file')
  setTogArgs(parser_tog)

  addCommonArgs(parser_split)
  addCommonArgs(parser_tog)

  return parser

def getRegexFiles(regex, from_fn):
  # Can specify a single file, in which case we should read it and return the
  # list of everything in the file.
  if from_fn:
    with open(regex) as f:
      return f.read().splitlines()

  import glob
  return glob.glob(regex)

def CleanName(prot_name):
  return prot_name.replace('.', '_').replace('/', '_')

# Dirty hack because multiprocessing can't call class functions.
def RunSingle(obj):
  return obj[0].runSingle(obj[1:])

# Class with functionality of running separate molecules.
class SeparateMols:
  def __init__(self, args):
    # Load the proteins.
    cmd.load(args.gold_r, 'pR')
    cmd.load(args.test_r, 'pRp')
    cmd.load(args.gold_l, 'pL')
    cmd.load(args.test_l, 'pLp')

    if args.prot_confs:
      self._lig_confs = getRegexFiles(args.prot_confs, args.input_fn)
      self._rec_confs = [_REC_TOGETHER] * len(self._lig_confs)
    else:
      self._lig_confs = getRegexFiles(args.lig_confs, args.input_fn)
      if args.rec_confs:
        self._rec_confs = getRegexFiles(args.rec_confs, args.input_fn)
      else:
        self._rec_confs = ['pRp.pdb'] * len(self._lig_confs)

  def GetLigandModels(self):
    return self._lig_confs

  def GetReceptorModels(self):
    return self._rec_confs

  def runSingle(self, x):
    lig = x[1]
    lig_name = CleanName(lig[:-4]) # Remove the '.pdb'
    cmd.load(lig, lig_name)


    # One case where these come together (SwarmDock).
    if x[2] == _REC_TOGETHER:
      return (lig_name, x[0].RMSDTog(lig_name))

    # Normal, separate function.
    rec = x[2]
    rec_name = CleanName(rec[:-4])
    if rec_name != 'pRp':
      cmd.load(rec, rec_name)

    rms = x[0].RMSDSep(rec_name, lig_name)
    #print 'Finished, rms is %f' % rms
    return (lig_name, rms)

class TogetherMols:
  def __init__(self, args):
    cmd.load(args.gold_p, 'goldp')
    cmd.load(args.test_p, 'testp')
    cmd.create('pR', 'goldp & chain %s' % args.gold_rec_chains)
    cmd.create('pL', 'goldp & chain %s' % args.gold_lig_chains)
    cmd.create('pRp', 'testp & chain %s' % args.test_rec_chains)
    cmd.create('pLp', 'testp & chain %s' % args.test_lig_chains)
    assert cmd.count_atoms('pR') > 0
    assert cmd.count_atoms('pL') > 0
    assert cmd.count_atoms('pRp') > 0
    assert cmd.count_atoms('pLp') > 0
    for obj in cmd.get_names():
      print 'obj is', obj
    self._gold_rec_chains = args.gold_rec_chains
    self._gold_lig_chains = args.gold_lig_chains
    self._test_rec_chains = args.test_rec_chains
    self._test_lig_chains = args.test_lig_chains

    self._test_confs = getRegexFiles(args.test_confs, args.input_fn)

  def GetLigandModels(self):
    return self._test_confs

  def GetReceptorModels(self):
    return [None]*len(self._test_confs)

  def runSingle(self, x):
    prot = x[1] # x[2] is None
    prot_name = CleanName(prot[:-4])
    cmd.load(prot, prot_name)

    rec_name = prot_name + '_R'
    lig_name = prot_name + '_L'

    cmd.select('%s & chain %s' %(prot_name, self._test_rec_chains), rec_name)
    cmd.select('%s & chain %s' %(prot_name, self._test_lig_chains), lig_name)

    rms = x[0].RMSDSep(lig_name, rec_name)
    return (prot_name, rms)

def RunMols(args, separated):
  if separated:
    runner = SeparateMols(args)
  else:
    runner = TogetherMols(args)

  # Need to do some cleanup of these.
  cmd.alter('all', 'segi=""')
  cmd.remove('not (alt ""+A)')
  cmd.alter('all', 'alt=""')

  alnr = PyMolAligner('pR', 'pRp', 'pL', 'pLp', args.dist)

  lig_files = runner.GetLigandModels()
  rec_files = runner.GetReceptorModels()
  if _VERBOSE:
    print 'ligs are', lig_files
    print 'recs are', rec_files
  assert len(lig_files) == len(rec_files)

  if not _VERBOSE:
    # Get rid of unwanted output.
    cmd.feedback('disable', 'executive', 'results')

  p = Pool(args.nproc)
  all_rms = p.map(RunSingle, [(runner, alnr, lig_files[i], rec_files[i]) for i in
    range(len(lig_files))])

  for i in range(len(lig_files)):
    print '%s %f' % all_rms[i]

def main():
  args = get_parser().parse_args()
  RunMols(args, hasattr(args, 'gold_r'))

if __name__ == 'pymol':
  main()
