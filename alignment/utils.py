""" Some utils to use other places. """
from __future__ import print_function

from Bio.SeqUtils import seq1
import os, sys, inspect
from pymol import cmd

_VERBOSE=0

def eassert(*asrt, **kwargs):
    line_num = ''.join([str(inspect.stack()[1][1]),
          ":", str(inspect.stack()[1][2]), ":", str(inspect.stack()[1][3])])
    fail_on_error = kwargs.get('fail_on_error', True)
    if not asrt[0]:
      msg = 'failure at line [' + str(line_num) + ']'
      if len(asrt) > 1:
        msg += ' : ' + asrt[1]
        if fail_on_error:
          raise AssertionError(msg)
        else:
          return False
    return True

def eprint(*args, **kwargs):
  if _VERBOSE:
    print(*args, file=sys.stderr, **kwargs)

def _get_unique_chain(chains, default_chain='A'):
  while default_chain in chains:
    default_chain = chr(ord(default_chain) + 1)

  return default_chain

def load_pdb(pdb_fn, pdb_str, default_chain='A'):
  # Delete it if it already exists to avoid appending state.
  cmd.delete(pdb_str)
  eprint('load %s, %s' %(pdb_fn, pdb_str))
  cmd.load(pdb_fn, pdb_str)

  # Need to do some cleanup of these.
  cmd.alter(pdb_str, 'segi=""')
  cmd.remove('%s & not (alt ""+A)' % pdb_str)
  cmd.alter(pdb_str, 'alt=""')

  # If there is no chain, we're going to break later. Just add something
  # default.
  chains = cmd.get_chains(pdb_str)
  eprint("%s: Chains are %s %d" % (pdb_str, chains, len(chains)))

  if '' in chains:
    next_chain = _get_unique_chain(chains, default_chain)
    # If there's only one chain, don't want this string.
    not_chain = ''
    if len(chains) > 1:
      # filter(None, chains) removes empty string.
      not_chain = '& not chain %s' % '+'.join(filter(None, chains))

    cmd.alter('%s %s' %(pdb_str, not_chain), 'chain="%s"' % next_chain)
    eprint("  after fix, chains are %s" % cmd.get_chains(pdb_str))

def getSequence(selection):
  myspace = {'residues': []}
  cmd.iterate('%s & n. ca' % selection, 'residues.append(resn)', space=myspace)
  return seq1(''.join(myspace['residues']))

def SimplifySelection(prot, cont):
  # Multiple chains requires two parentheses.
  multi_chain = False
  if len(set([x[0] for x in cont])) > 1:
    multi_chain = True

  str = '%s & n. ca & ' % prot
  if multi_chain:
    str += '('
  chain = None
  prev_resi = -2
  printed_resi = False
  for x in cont:
    this_chain = x[0]
    this_resi = int(x[2])
    # Initial chain not set.
    if not chain:
      str += '(chain %s & resi ' % this_chain
    # Found a new chain to process.
    elif chain != this_chain:
      if not printed_resi:
        str += '-%s' % prev_resi
      str += ') + (chain %s & resi ' % this_chain
      prev_resi = -2
      printed_resi = False

    # Found something non-sequential.
    if this_resi > prev_resi + 1:
      # Actually, haven't seen anything before (in this chain).
      if prev_resi == -2:
        str += '%s' % this_resi
      elif not printed_resi:
        # Haven't printed the terminal residue in this sequence.
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

  if multi_chain:
    str += ')'
  return str + ')'

