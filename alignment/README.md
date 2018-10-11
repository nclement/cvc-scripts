# Installation

This directory contains a few scripts that run with Pymol. One option is to run
these scripts with python, ensuring the Pymol library is in PYTHONPATH.
But...this doesn't work for the Educational Pymol. So we need to run it directly
with Pymol.... Which means we need to add a hack to include this directory in
the Pymol's Python path.

If you export this directory and add it to PYTHONPATH, it should work.

# Scripts Contained Herein

 * `alignRMSD_pymol_samechain.py`: Computes the cRMSD for a set of proteins,
   assuming they all come from the same chain (hack to make it work). See header
   of the script.
 * `alignRMSD_separate.py`: For a single protein (gold + test), computes cRMSD
   for separate (R&L) as well as full RMSD.
 * `alignRMSD_single.py`: Computes the one-sided cRMSD (only cRMSD for relevant
   atoms) for one or more proteins.

