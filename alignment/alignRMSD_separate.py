""" Compute the RMSD between a pair of proteins, including contact RMSD."""

from __future__ import print_function

from PyMolAligner import PyMolAligner
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

  parser.add_argument('-r', '--test_r', dest='test_r', required=True,
                      help='Test Receptor chain (or molecule)')
  parser.add_argument('-l', '--test_l', dest='test_l', required=True,
                      help='Test Ligand chain (or molecule)')

  # Vreven et al 2015 (zlab5) suggests to use 10\AA as the limit
  parser.add_argument('-X', '--dist', dest='dist', default=10, type=int,
                      help='Max distance between contact atoms')

  parser.add_argument('-V', '--verbose', dest='verbose', default=False,
      help='Print quite a bit of output',
      action='store_true')

  return parser

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

  utils.load_pdb(args.gold_r, 'pR')
  utils.load_pdb(args.test_r, 'pRp')
  utils.load_pdb(args.gold_l, 'pL')
  utils.load_pdb(args.test_l, 'pLp')

  alnr = PyMolAligner('pR', 'pRp', 'pL', 'pLp', args.dist)

  # Print Contact RMSD
  cont = alnr.RMSDSep('pRp', 'pLp', also_separate=True)
  full = alnr.FullRMSDSep('pRp', 'pLp')
  print("RMS is cont:%f %f %f" % cont, "full:%f %f %f" % full)

if __name__ == 'pymol':
  main()
