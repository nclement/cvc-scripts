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
prot=sys.argv[1]  # Real protein
samp=sys.argv[2]  # Sampled protein
chainsR=sys.argv[3]  # Receptor chains
chainsL=sys.argv[4]  # Ligand chains
# Vreven et al 2015 (zlab5) suggests to use 10\AA as the limit
X=10
if len(sys.argv) > 5:
  X=int(sys.argv[5])
#align_X=10
#align_X=5

verbose=0
if len(sys.argv) > 6:
  verbose=1


printVerbose("Running with args: %s %s %s %s [%d]" %(prot, samp, chainsR, chainsL, X))

cmd.load(prot, 'gold')
cmd.load(samp, 'testp')
cmd.alter('all', 'segi=""')
#select contact_rec_10, ((gold & chain r) & (all within $align_X of (gold and chain l))) & name ca 
#select contact_lig_10, ((gold & chain l) & (all within $align_X of (gold and chain r))) & name ca
cmd.select('contact_rec', '((gold & chain %s) & (all within %d of (gold and chain %s))) & n. ca' % (
        chainsR, X, chainsL))
cmd.select('contact_lig', '((gold & chain %s) & (all within %d of (gold and chain %s))) & n. ca' % (
        chainsL, X, chainsR))
cmd.select('contact_all', 'contact_rec + contact_lig')
contact_len = len(cmd.get_model('contact_all').atom)
# Do the actual alignment with no outlier rejection.
aln = cmd.super('testp', 'gold & contact_all', cycles=0)
method = 'super'
if aln[1] != contact_len:
  aln = cmd.align('testp', 'gold & contact_all', cycles=0, quiet=1)
  method = 'align'
if aln[1] != contact_len:
  aln = cmd.cealign('testp', 'gold & contact_all')
  aln = (aln['RMSD'], aln['alignment_length'], 'NA', 'NA', aln['alignment_length'])
  method = 'cealign'
if aln[1] != contact_len:
  fit_rms = cmd.fit('testp & n. ca', 'gold & contact_all', matchmaker=2, cycles=0, object='fit_aln')
  fit_atoms = int(len(cmd.get_model('fit_aln').atom)) / 2
  aln = (fit_rms, fit_atoms, 'NA', 'NA', fit_atoms)
  method='fit'
# Failed!
if aln[1] != contact_len:
  rms = cmd.rms('testp', 'gold & contact_all', cycles=0)
  printVerbose("RMS     RMSD: %f" % rms)
  print "RMS = NA NA"
else:
  print "RMS =", aln[0], method
