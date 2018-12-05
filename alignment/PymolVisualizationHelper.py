""" Output text that can be used to display cRMSD atoms. """

from PyMolAligner import PyMolAligner
import utils

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
  parser.add_argument('-V', '--verbose', dest='verbose', default=False,
      help='Print quite a bit of output',
      action='store_true')

  return parser

def main():
  try:
    args = get_parser().parse_args()
  except SystemExit, ext:
    return

  # Print the output from this.
  utils._VERBOSE=1
  utils.load_pdb(args.gold_r, 'pR')
  utils.load_pdb(args.test_r, 'pRp')
  utils.load_pdb(args.gold_l, 'pL')
  utils.load_pdb(args.test_l, 'pLp')

  if args.verbose:
    print('Using verbose output')
    utils._VERBOSE=1
  else:
    utils._VERBOSE=0

  alnr = PyMolAligner('pR', 'pRp', 'pL', 'pLp', args.dist)

  print (alnr)

if __name__ == "pymol":
  main()
