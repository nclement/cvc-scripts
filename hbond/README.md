# Installation

This directory contains git submodules. Make sure you run

```
git clone --recursive https://github.com/nclement/cvc-scripts
```
from the previous directory.

`protein_prep` _used_ to be a submodule, but I made some changes. So it's
included as a directory, but credit goes [here](https://github.com/Acpharis/protein_prep).

Everything in this directory doesn't (shouldn't) need to be compiled. Just
make sure to set all the programs in the top directory's Makefile.def file so
the correct binaries can be found.

# Background
Most of the programs used in this directory are based off programs developed
by the Structural Bioinformatics Library at Boston University. Original github
repositories can be found at https://github.com/StructuralBioinformaticsLab .

# Update Information
For your information, the following actions were performed.

## Parameters file
To update the parameters file, do the following:

```
curl http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_aug14.tgz | tar xz
cp toppar/toph19.inp prms/pdbamino_new.rtf
cat toppar/param19.inp toppar/ace/acepar19.inp > prms/parm_new.prm
```

The file `toopar/param19.inp` contains the effective Voronoi atom volumes for
each atom, and the file `toppar/ace/acepar19.inp` contains the parameters for
torsion angles, etc. (When discussed externally, this is the "charmm19 
forcefield"; other forcefields can be used.)

It is unclear where the atom prm file came from, but it can be downloaded from
[here](https://github.com/StructuralBioinformaticsLab/sidechain_example/blob/master/prms/atoms.0.0.6.prm.ms.3cap%2B0.5ace.Hr0rec)

The parm.prm file contains expected parameters for proteins. There are a few 
tricks:

 - If you see an error, "error on pdb atom %d: positive vdw eps for atom type %s",
   it's likely because the parm.prm file doesn't have the atom type defined
   in the NONBONDED section. You can create wildcards in this file, something
   like O\*, which will apply to all oxygens.
### Parameter file usage
There are several functions that use the parameters files. These are:

 - `read_psf(PSF_FILE) : atoms, atom_names, bonds, angles, dihedrals, donors, acceptors`:

    The PSF file is basically a PDB file with extra information computed by charmm (such
    as bonds, hbond donors and acceptors, etc). This function reads the following sections:

    1. `NATOM`: marker that designates the section for all atoms; needs to be
       in the "old" format, which is:

      ```
      ATOM_NUM SEG_NAME RES_ID RES_NAME ATOM_NAME ATOM_TYPE_NUM CHARGE ...
      ```
    2. `NBOND`: all of the bonded atoms, each line which has four pairs of
      atom numbers per line.
    3. `NTHETA`: theta angles for residues; lists three triples of atoms per
      line
    4. `NPHI` and `NIMPHI`: dihedrals and impropers, two quadruples of atoms
      per line
    5. `NDON` and `NACC`: donors and acceptors, four pairs of atoms per line

    After this has been done, the function compares string names and converts
    to residue numbers.

 - `read_atom_type_name_num(RTF_FILE) : atom_num -> atom_name mapping`:

    reads lines of the format: `MASS <ATOM_NUM> <ATOM_NAME>` to produce
    specific number to name mapping.

 - `get_atnam_rtf(RTF_FILE, ATOM_TYPE_NUMS) : atom_names`:

    returns a list of atom names from the (specific) atom types from input.

 - `read_pa(PRM_FILE, ATOM_NAME, BONDS, ANGLES, DIHEDRALS) : bond parameters`:

    The PRM file has all of the parameters for individual protein atoms, including
    torsion and bond angle expected values, and volumes of individual atoms. It
    is necessary that the values in this file and the values of the computed PSF
    file correlate.

## Additional Scripts
This folder also should contain the following additional submodules:

 - https://github.com/Acpharis/protein_prep
 - https://github.com/Acpharis/create_psf


## Common errors from running `testHbond`
| Error                                     | (Potential) Solution  |
| ----------------------------------------- | --------- |
| `KeyError: 'CD1'` from `create_psf` | Some versions of RTF files expect CD instead of CD1. There is a line in `rotein_prep/prepare.py` that will convert all of them for you. Later versions of RTF files expect CD1, so this line has been removed. If it causes you problems, try and fix by searching for `cd1_to_cd` in `protein_prep/prepare.py`. |
| `Warning: Ace Volume not set for atom %d` | It is likely that there is a disagreement between the `pdbamino.rtf` file and the `parm.prm` file. Make sure that for every atom name in the `pdbamino.rtf` file there is an associated volume in the VOLUME section of the `parm.prm` file. Alternatively, you're likely using a different `pdbamino.rtf` file for `protein_prep/prepare.py` and `testHbond`. |
| `error on pdb atom %d: positive vdw eps for atom type %s` | It's likely the `parm.prm` file doesn't have the atom type defined in the NONBONDED section. You can create wildcards for atoms (they should appear before all the non-wildcarded atoms), and should be of the form `O*` (applies to all Oxygens) and `O?` (applies to `OD` but not `OD1`) |

