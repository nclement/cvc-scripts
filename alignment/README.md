This directory contains a few scripts that run with Pymol. One option is to run
these scripts with python, ensuring the Pymol library is in PYTHONPATH.
But...this doesn't work for the Educational Pymol. So we need to run it directly
with Pymol.... Which means we need to add a hack to include this directory in
the Pymol's Python path.

If you export this directory and add it to PYTHONPATH, it should work.
