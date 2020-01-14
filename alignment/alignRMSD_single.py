""" Compute the contact RMSD of a single protein to its native structure

First, compute the contact RMSD between R+L, let this be {C_R+C_L}. Then, for
each R', compute the RMSD of R'*C_R and R*C_R. This is for seeing which sampled
protein is closer to the truth.

Will work with one or more of R' and L'.
"""

from __future__ import print_function

from PyMolAligner import PyMolAligner
from multiprocessing import Pool
import utils
from utils import eprint,eassert

def get_parser():
  """Set up arguments.
  
  Returns argparse object, e.g.:
    args = get_parser().parse_args()
    print args.dist
  """
  from argparse import ArgumentParser
  
  parser = ArgumentParser(description='Contact RMSD')

  parser.add_argument('-R', '--gold_r', dest='gold_r', required=True,
                      help='Gold Receptor chain (or molecule)')
  parser.add_argument('-L', '--gold_l', dest='gold_l', required=True,
                      help='Gold Ligand chain (or molecule)')

  parser.add_argument('-r', '--test_r', dest='test_r',
                      help='Test Receptor chain (or molecule)')
  parser.add_argument('-l', '--test_l', dest='test_l',
                      help='Test Ligand chain (or molecule)')
  parser.add_argument('--test_r_files', dest='test_r_files',
                      help='File with test Receptor molecules')
  parser.add_argument('--test_l_files', dest='test_l_files',
                      help='File with test Ligand molecules')

  # Vreven et al 2015 (zlab5) suggests to use 10\AA as the limit
  parser.add_argument('-X', '--dist', dest='dist', default=10, type=int,
                      help='Max distance between contact atoms')
  parser.add_argument('-N', dest='nproc', default=1, type=int,
      help='Number of processor cores to use')

  parser.add_argument('-V', '--verbose', dest='verbose', default=False,
      help='Print quite a bit of output',
      action='store_true')

  return parser

def get_test_structures(pp, p_list):
  if pp is None and p_list is None:
    return []

  if p_list is None:
    return [pp]

  prot_list = []
  with open(p_list) as f:
    prot_list = f.read().splitlines()
  if pp is not None:
    prot_list.append(pp)

  return prot_list

def CleanName(prot_name):
  return prot_name.replace('.', '_').replace('/', '_')

# Dirty hack because multiprocessing can't call class functions.
def RunSingleL(obj):
  lig = obj[2]
  lig_name = CleanName(lig[:-4])  # Remove the '.pdb'
  utils.load_pdb(lig, lig_name)
  return obj[0].RMSDSep(obj[1], lig_name, also_separate=True)
def RunSingleR(obj):
  rec = obj[1]
  rec_name = CleanName(rec[:-4])  # Remove the '.pdb'
  utils.load_pdb(rec, rec_name)
  return obj[0].RMSDSep(rec_name, obj[2], also_separate=True)

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

  # Some logic to prevent breaking when a p*p is not defined.
  pR = args.gold_r
  pL = args.gold_l
  pRp = args.gold_r
  pLp = args.gold_l
  pRList = get_test_structures(args.test_r, args.test_r_files)
  pLList = get_test_structures(args.test_l, args.test_l_files)
  if len(pRList) != 0:
    pRp = pRList[0]
  if len(pLList) != 0:
    pLp = pLList[0]

  utils.load_pdb(pR, 'pR')
  utils.load_pdb(pRp, 'pRp')
  utils.load_pdb(pL, 'pL')
  utils.load_pdb(pLp, 'pLp')

  alnr = PyMolAligner('pR', 'pRp', 'pL', 'pLp', args.dist)

  # Test with a single one because Pool hides traceback.
  if utils._VERBOSE:
    if len(pRList) != 0:
      rms = RunSingleR((alnr, pRList[0], 'pLp'))
      prot = pRList[0]
    else:
      rms = RunSingleL((alnr, 'pRp', pLList[0]))
      prot = pLList[0]
    print("%s %s" % (prot, rms))

  p = Pool(args.nproc)
  all_rms_l = p.map(RunSingleL, [(alnr, 'pRp', pLList[i]) for i in
    range(len(pLList))])
  all_rms_r = p.map(RunSingleR, [(alnr, pRList[i], 'pLp') for i in
    range(len(pRList))])

  # Print the RMSD here.
  for i, l in enumerate(pLList):
    print("%s %f" % (l, all_rms_l[i][2]))
  for i, r in enumerate(pRList):
    print("%s %f" % (r, all_rms_r[i][1]))

if __name__ == 'pymol':
  main()
