#! /usr/bin/perl
# This script is intented to fix the pdbqt file produced by openbabel
# to play more nicely with autodock vina. It's expecting a space to occur
# between numbers in the BRANCH and ENDBRANCH fields.
#
# I think this is the fix it's looking for.

while(<>) {
  if (/^BRANCH(....)(....)/) {
    print("BRANCH $1 $2\n");
  } elsif (/^ENDBRANCH(....)(....)/) {
    print("ENDBRANCH $1 $2\n");
  } else {
    print($_);
  }
}

