""" Make this a class so we can use it elsewhere. """
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from math import log
import pymol
import utils
from utils import eprint, eassert

def _getAlignmentList(str1, str2):
  """Get the alignment mapping from str1 (gold) => str2 (test). At the end,
  len(alns) == len(goldp)
  """
  alns = []
  idx1 = 0
  idx2 = 0
  for i in range(len(str1)):
    c1 = str1[i]
    c2 = str2[i]

    if c2 == '-':
      # Aligned with a gap in test; place a -1
      alns.append(-1)
      idx1 += 1
    elif c1 == '-':
      # Aligned with a gap in gold; just skip.
      idx2 += 1
      continue
    else:
      # Good on both top and bottom.
      alns.append(idx2)
      idx1 += 1
      idx2 += 1
  return alns

def _testContactRes(aln_list, contact_idx, full_seq):
  """From the alignment list, aln_list, get the contact residue indices for
  the test seq (similar to getContactResIndx for gold).

  Arguments:
    aln_list: alignment list from getAlignmentList
    contact_idx: contact residue indices for gold
    full_seq: the tuple (from getFullRes) representing the full sequence
       for this test protein.
  
  """
  test_cont = []
  for x in contact_idx:
    # This means it aligns to nothing--will cause the number of atoms in
    # test and gold to be different (bad for RMS).
    eassert(aln_list[x] != -1,
        'shouldn\'t align to nothing! at pos %d:%d' %(x, aln_list[x]))

    #eprint(' - test appending x=%d, len %d,%d'
    #    %(x, len(full_seq), len(aln_list)))
    #eprint('   which is', aln_list[x], full_seq[aln_list[x]])
    test_cont.append(full_seq[aln_list[x]])
  return test_cont

def _gap_function(x, y):  # x is gap position in seq, y is gap length
  if y == 0:  # No gap
    return 0
  elif y == 1:  # Gap open penalty
      return -2
  return - (2 + y/4.0 + log(y)/2.0)
blosum_matrix = matlist.blosum62
def _getBestAlignment(seq1, seq2):
  eprint('Doing the alignment...')
  alns = pairwise2.align.globaldc(seq1, seq2, blosum_matrix,
      _gap_function, _gap_function, one_alignment_only=True)
  eprint('Finished!')
  return alns[0]

def _getContactRes(sel1, sel2, X):
  pymol.cmd.select('contact_rec', '(%s) & (all within %s of (%s)) & n. ca'
      %(sel1, X, sel2))
  sp = {'contact': []}
  pymol.cmd.iterate('contact_rec', 'contact.append((chain, resn, resi, index))', space=sp)
  return sp['contact']

def _getFullRes(sel):
  sp = {'full': []}
  pymol.cmd.iterate('%s & n. ca' % sel, 'full.append((chain, resn, resi, index))',
      space=sp)
  return sp['full']

def _fixContact(contact, contact_idx, aln_list, fatal_missing):
  # Remove things that align with nothing in the other sequence.
  no_align = []
  for i, x in enumerate(contact_idx):
    if aln_list[x] == -1:
      no_align.append(i)
  for i, x in enumerate(no_align):
    # Each thing we delete will reduce the position by 1.
    eprint('Deleting contact at position %d=%d' % (x, contact_idx[x-i]))
    if fatal_missing:
      eassert(False, 'Error aligning to a gap! Remove --fail_on_missing_atoms '
          'to continue.')
      
    del contact[x-i]
    del contact_idx[x-i]

  return (contact, contact_idx)


def _getContactResIndx(contact, full):
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


class PyMolAligner:
  def __init__(self, protR, protRp, protL, protLp, X,
      fail_on_missing_atoms=False):
    # Remember these.
    self._protR = protR
    self._protL = protL
    self._protRp = protRp
    self._protLp = protLp

    self._gold_chains_r = pymol.cmd.get_chains(protR)
    self._gold_chains_l = pymol.cmd.get_chains(protL)
    self._test_chains_r = pymol.cmd.get_chains(protRp)
    self._test_chains_l = pymol.cmd.get_chains(protLp)

    # Get the sequences.
    self._seqR_gold = utils.getSequence(protR)
    self._seqL_gold = utils.getSequence(protL)
    self._seqR_test = utils.getSequence(protRp)
    self._seqL_test = utils.getSequence(protLp)

    eprint('R_gold is', self._seqR_gold)
    eprint('R_test is', self._seqR_test)
    self._R_aln = _getBestAlignment(self._seqR_gold, self._seqR_test)
    self._L_aln = _getBestAlignment(self._seqL_gold, self._seqL_test)
    eprint('R_aln is', self._R_aln)
    eprint('L_aln is', self._L_aln)

    self._full_R = _getFullRes(protR)
    self._full_L = _getFullRes(protL)
    self._full_Rp = _getFullRes(protRp)
    self._full_Lp = _getFullRes(protLp)
    
    self._aln_list_R = _getAlignmentList(self._R_aln[0], self._R_aln[1])
    self._aln_list_L = _getAlignmentList(self._L_aln[0], self._L_aln[1])

    # Need to get contact residues for gold, then corresponding res from test.
    self._cont_R = _getContactRes(protR, protL, X)
    self._cont_L = _getContactRes(protL, protR, X)

    self._cont_R_idx = _getContactResIndx(self._cont_R, self._full_R)
    self._cont_L_idx = _getContactResIndx(self._cont_L, self._full_L)

    eprint('Before fixing, R contact is', self._cont_R, self._cont_R_idx)
    eprint('Before fixing, L contact is', self._cont_L, self._cont_L_idx)
    # If there are gaps in the alignment, should skip these (otherwise the
    # RMS calculation will be incorrect).
    (self._cont_R, self._cont_R_idx) = _fixContact(
        self._cont_R, self._cont_R_idx, self._aln_list_R, fail_on_missing_atoms)
    (self._cont_L, self._cont_L_idx) = _fixContact(
        self._cont_L, self._cont_L_idx, self._aln_list_L, fail_on_missing_atoms)

    eprint('R Contact is', self._cont_R)
    eprint('R Contact_idx is', self._cont_R_idx)
    eprint('R aln is', self._aln_list_R)
    eprint('R full is', self._full_Rp)
    eprint('L Contact is', self._cont_L)
    eprint('L Contact_idx is', self._cont_L_idx)
    eprint('L aln is', self._aln_list_L)
    eprint('L full is', self._full_Lp)
    self._cont_Rp = _testContactRes(self._aln_list_R, self._cont_R_idx, self._full_Rp)
    self._cont_Lp = _testContactRes(self._aln_list_L, self._cont_L_idx, self._full_Lp)
    eprint('R cont (len:%d) is' % len(self._cont_R), self._cont_R)
    eprint('L cont (len:%d) is' % len(self._cont_L), self._cont_L)
    eprint('Rp cont (len:%d) is' % len(self._cont_Rp), self._cont_Rp)
    eprint('Lp cont (len:%d) is' % len(self._cont_Lp), self._cont_Lp)

  def __str__(self):
    retstr = ''
    # Load objects.
    #retstr += 'load %s as gold_r\n' % self._protR
    #retstr += 'load %s as gold_l\n' % self._protL
    #retstr += 'load %s as test_r\n' % self._protRp
    #retstr += 'load %s as test_l\n' % self._protLp
    # Select contact atoms.
    simp_R = utils.SimplifySelection(self._protR, self._cont_R)
    simp_Rp = utils.SimplifySelection(self._protRp, self._cont_Rp)
    simp_L = utils.SimplifySelection(self._protL, self._cont_L)
    simp_Lp = utils.SimplifySelection(self._protLp, self._cont_Lp)
    retstr += 'select cont_gold_r, %s\n' % simp_R
    retstr += 'select cont_gold_l, %s\n' % simp_L
    retstr += 'select cont_test_r, %s\n' % simp_Rp
    retstr += 'select cont_test_l, %s\n' % simp_Lp
    retstr += 'align cont_test_r, cont_gold_r\n'
    retstr += 'align cont_test_l, cont_gold_l\n'
    retstr += 'show_as cartoon\n'
    return retstr

  def RMSDTog(self, prot):
    """Returns the RMSD between gold and test, specifying a single structure and
    getting R+L from chains."""
    simp_R = utils.SimplifySelection(self._protR, self._cont_R)
    simp_Rp = utils.SimplifySelection(prot, self._cont_Rp)
    simp_L = utils.SimplifySelection(self._protL, self._cont_L)
    simp_Lp = utils.SimplifySelection(prot, self._cont_Lp)
    eprint('%s: %d vs %d' % (simp_R, pymol.cmd.count_atoms(simp_R),
      len(self._cont_R)))
    eprint('%s: %d vs %d' % (simp_Rp, pymol.cmd.count_atoms(simp_Rp),
        len(self._cont_Rp)))
    eprint('%s: %d vs %d' % (simp_L, pymol.cmd.count_atoms(simp_L),
      len(self._cont_L)))
    eprint('%s: %d vs %d' % (simp_Lp, pymol.cmd.count_atoms(simp_Lp),
        len(self._cont_Lp)))
    eassert(pymol.cmd.count_atoms(simp_R) == len(self._cont_R),
        'number of atoms in contact selection bad for gold:R')
    eassert(pymol.cmd.count_atoms(simp_Rp) == len(self._cont_Rp),
        'number of atoms in contact selection bad for test:R')
    eassert(pymol.cmd.count_atoms(simp_L) == len(self._cont_L),
        'number of atoms in contact selection bad for gold:L')
    eassert(pymol.cmd.count_atoms(simp_Lp) == len(self._cont_Lp),
        'number of atoms in contact selection bad for test:L')
    return pymol.cmd.pair_fit(simp_R, simp_Rp, simp_L, simp_Lp)

  def RMSDSep(self, rec, lig, also_separate=False, fail_on_error=True):
    """Returns the RMSD between two proteins when they are specified as two
    objects.

    Passing also_separate=True will return a tuple, which is:
      (contact_rmsd_RL, contact_rmsd_R, contact_RMSD_L)
    """
    simp_R = utils.SimplifySelection(self._protR, self._cont_R)
    simp_Rp = utils.SimplifySelection(rec, self._cont_Rp)
    simp_L = utils.SimplifySelection(self._protL, self._cont_L)
    simp_Lp = utils.SimplifySelection(lig, self._cont_Lp)
    eprint('%s: %d vs %d' % (simp_R, pymol.cmd.count_atoms(simp_R),
      len(self._cont_R)))
    eprint('%s: %d vs %d' % (simp_Rp, pymol.cmd.count_atoms(simp_Rp),
      len(self._cont_Rp)))
    eprint('%s: %d vs %d' % (simp_L, pymol.cmd.count_atoms(simp_L),
      len(self._cont_L)))
    eprint('%s: %d vs %d' % (simp_Lp, pymol.cmd.count_atoms(simp_Lp),
      len(self._cont_Lp)))
    if not  eassert(pymol.cmd.count_atoms(simp_R) == len(self._cont_R),
        'number of atoms in contact selection bad for gold:R', fail_on_error=fail_on_error):
      return None
    if not eassert(pymol.cmd.count_atoms(simp_Rp) == len(self._cont_Rp),
        'number of atoms in contact selection bad for test:R', fail_on_error=fail_on_error):
      return None
    if not eassert(pymol.cmd.count_atoms(simp_L) == len(self._cont_L),
        'number of atoms in contact selection bad for gold:L', fail_on_error=fail_on_error):
      return None
    if not eassert(pymol.cmd.count_atoms(simp_Lp) == len(self._cont_Lp),
        'number of atoms in contact selection bad for test:L', fail_on_error=fail_on_error):
      return None

    if also_separate:
      return (pymol.cmd.pair_fit(simp_R, simp_Rp, simp_L, simp_Lp),
          pymol.cmd.pair_fit(simp_R, simp_Rp),
          pymol.cmd.pair_fit(simp_L, simp_Lp))
    return pymol.cmd.pair_fit(simp_R, simp_Rp, simp_L, simp_Lp)

  # Compute the RMSD between all atoms, and not just the contact ones.
  def FullRMSDSep(self, rec, lig):
    full_idx_R = []
    full_sel_R = []
    for i, x in enumerate(self._aln_list_R):
      if x != -1:
        full_idx_R.append(i)
        full_sel_R.append(self._full_R[i])

    full_idx_L = []
    full_sel_L = []
    for i, x in enumerate(self._aln_list_L):
      if x != -1:
        full_idx_L.append(i)
        full_sel_L.append(self._full_L[i])

    full_sel_Rp = _testContactRes(self._aln_list_R, full_idx_R, self._full_Rp)
    full_sel_Lp = _testContactRes(self._aln_list_L, full_idx_L, self._full_Lp)

    simp_R = utils.SimplifySelection(self._protR, full_sel_R)
    simp_Rp = utils.SimplifySelection(rec, full_sel_Rp)
    simp_L = utils.SimplifySelection(self._protL, full_sel_L)
    simp_Lp = utils.SimplifySelection(lig, full_sel_Lp)
    eprint('%s: %d vs %d' % (simp_R, pymol.cmd.count_atoms(simp_R),
      len(full_sel_R)))
    eprint('%s: %d vs %d' % (simp_Rp, pymol.cmd.count_atoms(simp_Rp),
      len(full_sel_Rp)))
    eprint('%s: %d vs %d' % (simp_L, pymol.cmd.count_atoms(simp_L),
      len(full_sel_L)))
    eprint('%s: %d vs %d' % (simp_Lp, pymol.cmd.count_atoms(simp_Lp),
      len(full_sel_Lp)))
    eassert(pymol.cmd.count_atoms(simp_R) == len(full_sel_R),
        'number of atoms in contact selection bad for gold:R')
    eassert(pymol.cmd.count_atoms(simp_Rp) == len(full_sel_Rp),
        'number of atoms in contact selection bad for test:R')
    eassert(pymol.cmd.count_atoms(simp_L) == len(full_sel_L),
        'number of atoms in contact selection bad for gold:L')
    eassert(pymol.cmd.count_atoms(simp_Lp) == len(full_sel_Lp),
        'number of atoms in contact selection bad for test:L')

    return (pymol.cmd.pair_fit(simp_R, simp_Rp, simp_L, simp_Lp),
        pymol.cmd.pair_fit(simp_R, simp_Rp),
        pymol.cmd.pair_fit(simp_L, simp_Lp))
