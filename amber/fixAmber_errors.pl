#! /usr/bin/perl

#
# Sometimes, tleap isn't able to generate the output because there are 
# unidentified hydrogen atoms (named wrong, etc). One solution is just to
# delete these atoms from the PDB and see if that works. We'll try it.
#

use strict;

my $pdb = shift;
my $t_error = shift;
my $pdb_out = shift; # Can be same as out.

die "usage: $0 <pdb> <tleap_error> <pdb_out>\n" if !$pdb || !$t_error || !$pdb_out;

open (my $PDB_IN, "<", $pdb) or die "Error: Could not open PDB input file: $!\n";
open (my $T_ERR, "<", $t_error) or die "Error: Could not open tleap error file: $!\n";

my @pdb_lines = <$PDB_IN>;
my @t_lines = <$T_ERR>;
my @errs = grep /FATAL:\s+Atom/, @t_lines;

# Close the files
close($PDB_IN);
close($T_ERR);

my %fixed = ();

for (@errs) {
  /.R<.?(...) (\d+)>.A<(\S+)/;
  my $res = $1;
  my $resi = $2;
  my $atom = $3;

  my $found_err = 0;
  for (my $i = 0; $i < @pdb_lines; ++$i) {
    # If it matches,
    # Ignore the position
    if ($pdb_lines[$i] =~ /$atom\s+$res/) {
      # Delete it
      $pdb_lines[$i] = "";
      $found_err = 1;
      #last; # Found, finished.
    }
  }
  # It's possible the residue is named something funky. So let's try and just
  # delete all the offending atoms.
  if (!$found_err) {
    # Ignore it if we've already fixed it.
    if (defined $fixed{$atom}) {
      print "Already fixed atom $atom, continuing\n";
      next;
    }

    print "Uh-oh. Couldn't remove error. Removing all matches to $atom\n";
    my $del = 0;
    for (my $i = 0; $i < @pdb_lines; ++$i) {
      my $this_atom = substr ($pdb_lines[$i], 13, 4);
      $this_atom =~ s/\s*(\S+)\s*/$1/;
      if ($this_atom eq $atom) {
        $pdb_lines[$i] = "";
        $del++;
        $found_err = 1;
        $fixed{$atom} = 1;
      }
    }
    print "  deleted $del things\n";
  }
  # If we still haven't found the error, it's possible that the atom name got
  # converted from something like HB1 to 1HB. See if we can strip off the
  # digits and try again.
  if (!$found_err) {
    $atom =~ s/^\d*([^\d]+)\d*$/$1/;
    # We've already dealt with this atom.
    if (defined $fixed{$atom}) {
      print "Already fixed atom $atom, continuing\n";
      next;
    }
    print "Hmm... Still can't find error. trying with atom $atom\n";
    my $del = 0;
    for (my $i = 0; $i < @pdb_lines; ++$i) {
      my $this_atom = substr ($pdb_lines[$i], 13, 4);
      $this_atom =~ s/^\s*\d*([^\d\s]+)\d*\s*$/$1/;
      #print "  [$atom] vs [$this_atom]\n";
      if ($this_atom eq $atom) {
        $pdb_lines[$i] = "";
        $del++;
        $found_err = 1;
        $fixed{$atom} = 1;
      }
    }
    print "  deleted $del things\n";
  }
  # Ooops.
  if (!$found_err) {
    print "Couldn't solve the error. Sorry.\n";
  }
}

open (my $PDB_OUT, ">", $pdb_out) or die "Error: could not open outfile: $!\n";
for (@pdb_lines) {
  next if /^$/;
  print $PDB_OUT $_;
}

