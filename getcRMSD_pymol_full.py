import sys
import os

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

verbose=1
if len(sys.argv) > 6:
  verbose=1

def renumber(pR, pL):
  num_rec = cmd.count_atoms('%s & n. ca' % pR)
  cmd.alter(pL, 'resi=str(int(resi)+%d)' % num_rec)

printVerbose("Running with args: %s %s %s %s [%d]" %(protR, protRp, protL, protLp, X))

cmd.load(protR, 'pR')
cmd.load(protRp, 'pRp')
cmd.load(protL, 'pL')
cmd.load(protLp, 'pLp')
cmd.alter('all', 'segi=""')
gold_chains_r = cmd.get_chains('pR')
gold_chains_l = cmd.get_chains('pL')
#test_chains_r = cmd.get_chains('pRp')
#test_chains_l = cmd.get_chains('pLp')

# I think we need to renumber these first...
renumber('pR', 'pL')
renumber('pRp', 'pLp')
cmd.create('gold', 'pL + pR')
cmd.create('testp', 'pLp + pRp')
#select contact_rec_10, ((gold & chain r) & (all within $align_X of (gold and chain l))) & name ca 
#select contact_lig_10, ((gold & chain l) & (all within $align_X of (gold and chain r))) & name ca
cmd.select('contact_rec', '((gold & chain %s) & (all within %s of (gold and chain %s))) & n. ca' 
    % ('+'.join(gold_chains_r), X, '+'.join(gold_chains_l)))
cmd.select('contact_lig', '((gold & chain %s) & (all within %s of (gold and chain %s))) & n. ca' 
    % ('+'.join(gold_chains_l), X, '+'.join(gold_chains_r)))
#cmd.select('contact_lig', '((gold & chain %s) & (all within %s of (gold and chain %s))) & n. ca' % X)
cmd.select('contact_all', 'contact_rec + contact_lig')
printVerbose('contact_all & n.ca length is', cmd.count_atoms('contact_all & n. ca'))
printVerbose('gold & n. ca length is: ', len(cmd.get_model('gold & n. ca').atom))
printVerbose('testp & n. ca length is: ', len(cmd.get_model('testp & n. ca').atom))
printVerbose('contact_rec & n. ca length is: ', len(cmd.get_model('contact_rec & n. ca').atom))
printVerbose('contact_lig & n. ca length is: ', len(cmd.get_model('contact_lig & n. ca').atom))
printVerbose('length of contact_all is: ', len(cmd.get_model('contact_all & n. ca').atom))
contact_len = len(cmd.get_model('contact_all').atom)
# Do the actual alignment with no outlier rejection.
# First, try super. Seems to work well.
aln = cmd.super('testp', 'gold & contact_all', cycles=0, quiet=1)
method='super'
printVerbose("Super   RMSD: %f, aligned after: %d aligned before: %d" % (aln[0], aln[1], aln[4]))

# If that doesn't work, try align. Less optimal.
if aln[1] != contact_len:
  aln = cmd.align('testp', 'gold & contact_all', cycles=0, quiet=1)
  method='align'
  printVerbose("ALIGN   RMSD: %f, aligned after: %d aligned before: %d" % (aln[0], aln[1], aln[4]))
if aln[1] != contact_len:
  aln = cmd.cealign('testp', 'gold & contact_all')
  aln = (aln['RMSD'], aln['alignment_length'], 'NA', 'NA', aln['alignment_length'])
  method='cealign'
  printVerbose("CEALIGN RMSD: %f, aligned after: %d aligned before: %d" % (aln[0], aln[1], aln[4]))
if aln[1] != contact_len:
  fit_rms = cmd.fit('testp & n. ca', 'gold & contact_all', matchmaker=2, cycles=0, object='fit_aln')
  fit_atoms = int(len(cmd.get_model('fit_aln').atom)) / 2
  aln = (fit_rms, fit_atoms, 'NA', 'NA', fit_atoms)
  method='fit'
  printVerbose("FIT     RMSD: %f, aligned after: %d aligned before: %d" % (aln[0], aln[1], aln[4]))
# Failed!
if aln[1] != contact_len:
  rms = cmd.rms('testp', 'gold & contact_all', cycles=0)
  printVerbose("RMS     RMSD: %f" % rms)
  print "RMS = NA NA"
else:
  print "RMS =", aln[0], method
