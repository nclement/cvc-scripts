""" Output text that can be used to display cRMSD atoms. """

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

  return parser

def main():
  args = get_parser().parse_args()
  alnr = PyMolAligner(args.gold_r, args.test_r, args.gold_l, args.test_l,
      args.dist)

  print (alnr)

if __name__ == "pymol":
  main()
