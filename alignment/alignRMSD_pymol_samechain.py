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


There are a few known issues:
  1. The alignment treats the entire protein as a single chain (i.e. single
  sequence). If there are more than one chain, the alignment might be wrong.
  This might be an issue, but we'll ignore it for now.
  Some fixes:
    a. Perform an alignment between each chain, and use DP to find the
       optimal chain-chain alignment. Instead, we're expecting the chains
       to come in the same order, as this solution would require extensive
       rewriting (consider each chain separately).

  2. We're only considering the top alignment. It's possible this will cause
  contact residues in the gold sequence to align with gaps in the test
  sequence. What should we do in this scenario? It was decided that we'd just
  ignore these residues, which might not be the best answer. Some fixes:
    a. Consider more than one alignment, see if there is something where
       everything fits.
    b. Perform only a local alignment, which might remove some of the gaps at
       the end.
"""
from contextlib import contextmanager
from multiprocessing import Pool
from PyMolAligner import PyMolAligner
import utils
from utils import eprint,eassert

_REC_TOGETHER="__tog"


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
  parser.add_argument('--fail_on_missing_atoms', action='store_true',
      help='If a contact atom is aligned to a gap, fail and exit. If this flag '
           'is not specified, the program will silently continue upond finding '
           'this case.')
  parser.add_argument('-V', '--verbose', dest='verbose', default=False,
      help='Print quite a bit of output',
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
    utils.load_pdb(args.gold_r, 'pR')
    utils.load_pdb(args.test_r, 'pRp')
    utils.load_pdb(args.gold_l, 'pL')
    utils.load_pdb(args.test_l, 'pLp')

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
    utils.load_pdb(lig, lig_name)


    # One case where these come together (SwarmDock).
    if x[2] == _REC_TOGETHER:
      return (lig_name, x[0].RMSDTog(lig_name))

    # Normal, separate function.
    rec = x[2]
    rec_name = CleanName(rec[:-4])
    if rec_name != 'pRp':
      utils.load_pdb(rec, rec_name)

    rms = x[0].RMSDSep(rec_name, lig_name)
    #print 'Finished, rms is %f' % rms
    return (lig_name, rms)

class TogetherMols:
  def __init__(self, args):
    utils.load_pdb(args.gold_p, 'goldp')
    utils.load_pdb(args.test_p, 'testp')

    eprint('goldp & chain %s' % '+'.join(args.gold_rec_chains))
    cmd.create('pR', 'goldp & chain %s' % '+'.join(args.gold_rec_chains))
    cmd.create('pL', 'goldp & chain %s' % '+'.join(args.gold_lig_chains))
    cmd.create('pRp', 'testp & chain %s' % '+'.join(args.test_rec_chains))
    cmd.create('pLp', 'testp & chain %s' % '+'.join(args.test_lig_chains))
    eassert(cmd.count_atoms('pR') > 0, 'no atoms in gold:R')
    eassert(cmd.count_atoms('pL') > 0, 'no atoms in gold:L')
    eassert(cmd.count_atoms('pRp') > 0, 'no atoms in test:R')
    eassert(cmd.count_atoms('pLp') > 0, 'no atoms in test:L')
    for obj in cmd.get_names():
      eprint('obj is %s' % obj)
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
    utils.load_pdb(prot, prot_name)

    rec_name = prot_name + '_R'
    lig_name = prot_name + '_L'

    eprint('cmd.select %s, %s & chain %s' % (
      rec_name, prot_name, '+'.join(self._test_rec_chains)))
    eprint('cmd.select %s, %s & chain %s' % (
      lig_name, prot_name, '+'.join(self._test_lig_chains)))
    cmd.select(rec_name, '%s & chain %s' %(
      prot_name, '+'.join(self._test_rec_chains)))
    cmd.select(lig_name, '%s & chain %s' %(
      prot_name, '+'.join(self._test_lig_chains)))

    rms = x[0].RMSDSep(rec_name, lig_name)
    return (prot_name, rms)

def RunMols(args, separated):
  if separated:
    runner = SeparateMols(args)
  else:
    runner = TogetherMols(args)

  alnr = PyMolAligner('pR', 'pRp', 'pL', 'pLp', args.dist,
      args.fail_on_missing_atoms)

  lig_files = runner.GetLigandModels()
  rec_files = runner.GetReceptorModels()
  eprint('ligs are', lig_files)
  eprint('recs are', rec_files)
  eassert(len(lig_files) == len(rec_files),
      'you didn\'t specify the correct number of rec/lig files!')

  # Test with a single one because Pool hides traceback.
  if utils._VERBOSE:
    RunSingle((runner, alnr, lig_files[0], rec_files[0]))

  p = Pool(args.nproc)
  all_rms = p.map(RunSingle, [(runner, alnr, lig_files[i], rec_files[i]) for i in
    range(len(lig_files))])

  for i in range(len(lig_files)):
    print('%s %f' % all_rms[i])

def main():
  if not utils._VERBOSE:
    # Get rid of unwanted output.
    cmd.feedback('disable', 'executive', 'results')

  # Hack to return cleanly instead of calling sys.exit() so Pymol can print
  # decent information.
  try:
    args = get_parser().parse_args()
  except SystemExit, ext:
    return

  if args.verbose:
    print('Using verbose output')
    utils._VERBOSE=1
  RunMols(args, hasattr(args, 'gold_r'))

if __name__ == 'pymol':
  main()
