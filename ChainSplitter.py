# Modification of a file from
# https://stackoverflow.com/questions/11685716/how-to-extract-chains-from-a-pdb-file
from __future__ import print_function

import os
import re
from Bio import PDB
from Bio.PDB import PDBList


class ChainSplitter:
    def __init__(self, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        self.parser = PDB.PDBParser()
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.getcwd()
        self.out_dir = out_dir

    def make_pdb(self, pdb_fn, chain_letters, pdb_id, out_path, overwrite=True):
        """ Create a new PDB file containing only the specified chains.

        Returns the path to the created file.

        :param pdb: full path to the crystal structure
        :param chain_letters: iterable of chain characters (case insensitive)
        :param overwrite: write over the output file if it exists
        """
        # If chain looks like 'A(6)', it means 'chain A model 6'.
        model = None
        p = re.match(r'(.*)\((\d+)\)', ''.join(chain_letters))
        if p:
          chain_letters = list(p.group(1))
          model = int(p.group(2))

        chain_letters = [chain.upper() for chain in chain_letters]

        # Input/output files
        struct = self.parser.get_structure('struct', pdb_fn)
        plural = "s" if (len(chain_letters) > 1) else ""  # for printing

        # Skip PDB generation if the file already exists
        if (not overwrite) and (os.path.isfile(out_path)):
            print("Chain%s %s of '%s' already extracted to '%s'." %
                    (plural, ", ".join(chain_letters), pdb_id, out_name))
            return out_path

        print("Extracting chain%s %s from %s..." % (plural,
                ", ".join(chain_letters), pdb_id))

        for c in chain_letters:
          if c == '?':
            continue

          if c not in [ch.id for ch in struct.get_chains()]:
            print("Error: Could not find chain %s in pdb, pdb has %s" % (
              c, [ch.id for ch in struct.get_chains()]))
            return None

        # Write new file with only given chains
        self.writer.set_structure(struct)
        self.writer.save(out_path, select=SelectChains(chain_letters, model))

        # Then add header information.
        goodLines = []
        with open(pdb_fn, 'r') as header, open(out_path, 'r') as atoms:
          for line in header:
            if (line.startswith("HELIX") or line.startswith("SHEET")
                or line.startswith("REMARK")):
              goodLines.append(line)
          for line in atoms:
            goodLines.append(line)

        with open(out_path, 'w') as outf:
          for line in goodLines:
            outf.write(line)

        return out_path


class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters, model):
        self.chain_letters = chain_letters
        self.model = model

    def accept_chain(self, chain):
      if self.chain_letters == ['?']:
        return True
      return (chain.get_id() in self.chain_letters)

    def accept_atom(self, atom):
      return (atom.get_altloc() in [' ', 'A'])

    def accept_model(self, model):
      if self.model:
        return model.serial_num == self.model
      return True


if __name__ == "__main__":
    """ Parses PDB id's desired chains, and creates new PDB structures. """
    import sys
    if not len(sys.argv) == 4:
        print("Usage: $ python %s pdb chains outf" % __file__)
        sys.exit()

    splitter = ChainSplitter()

    pdb_id = sys.argv[1]
    chains = sys.argv[2]
    outf = sys.argv[3]

    pdbList = PDBList()
    pdb_fn = pdbList.retrieve_pdb_file(
        pdb_id, file_format="pdb", pdir="tmp_files")
    print("Saved original to %s" % pdb_fn)
    outf = splitter.make_pdb(pdb_fn, list(chains), pdb_id, outf)
    if outf:
      print("Saved to %s" % outf)
