import sys

# Given a protein pair, find the contact-residue RMSD (cRMSD), which is 
# defined as all residues within X\AA of the docked protein complex.
#
# This is done as follows:
#  - Find all residues of A within X\AA of B
#  - Align A with A' along contact RMSD atoms
#  - Compute RMSD over all Ca atoms

def printVerbose(*str):
  if verbose:
    print str

if len(sys.argv) < 4:
  print("usage getcRMSD_pymol_full.py <GOLD_R> <TEST_R> <GOLD_L> <TEST_L> [LIMIT=10]")
  os._exit(1)
protR=sys.argv[1]  # Receptor
protRp=sys.argv[2] # sampled receptor
protL=sys.argv[3]  # Ligand
protLp=sys.argv[4] # sampled ligand
# Vreven et al 2015 (zlab5) suggests to use 10\AA as the limit
X=10
if len(sys.argv) > 5:
  X=int(sys.argv[5])
#align_X=10
#align_X=5

verbose=0
if len(sys.argv) > 6:
  verbose=1


printVerbose("Running with args: %s %s %s %s [%d]" %(protR, protRp, protL, protLp, X))

cmd.load(protR, 'pR')
cmd.load(protRp, 'pRp')
cmd.load(protL, 'pL')
cmd.load(protLp, 'pLp')
cmd.alter('all', 'segi=""')
cmd.alter('pR', 'chain="r"')
cmd.alter('pRp', 'chain="r"')
cmd.alter('pL', 'chain="l"')
cmd.alter('pLp', 'chain="l"')
cmd.create('gold', 'pL + pR')
cmd.create('testp', 'pLp + pRp')
#select contact_rec_10, ((gold & chain r) & (all within $align_X of (gold and chain l))) & name ca 
#select contact_lig_10, ((gold & chain l) & (all within $align_X of (gold and chain r))) & name ca
cmd.select('contact_rec', '((gold & chain r) & (all within %s of (gold and chain l))) & n. ca' % X)
cmd.select('contact_lig', '((gold & chain l) & (all within %s of (gold and chain r))) & n. ca' % X)
cmd.select('contact_all', 'contact_rec + contact_lig')
printVerbose('other length is', cmd.count_atoms('contact_all & n. ca'))
printVerbose('length is: ', len(cmd.get_model('gold & n. ca').atom))
printVerbose('length is: ', len(cmd.get_model('testp & n. ca').atom))
printVerbose('length is: ', len(cmd.get_model('contact_rec & n. ca').atom))
printVerbose('length is: ', len(cmd.get_model('contact_lig & n. ca').atom))
printVerbose('length is: ', len(cmd.get_model('contact_all & n. ca').atom))
# Do the actual alignment with no outlier rejection.
aln = cmd.super('testp', 'gold & contact_all', cycles=0)
printVerbose("RMSD: %f, aligned after: %d aligned before: %d" % (aln[0], aln[1], aln[4]))
print aln[0]
