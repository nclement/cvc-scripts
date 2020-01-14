#! /usr/bin/perl

use strict;

###############################################################################
# Want to be able to impose a global ordering on a PDB file so we can compare #
# them directly as long as they have the same atoms.                          #
# 
# This is the protocol we'll use:
#   1. Sort by residue number
#   2. Sort by atom ID#, which is assigned as follows:
#      N : 0
#      CA: 1
#      C : 2
#      O : 3
#      CB: 4
#      CG: 5
#      CD: 6
#      everything else: 7
#   3. Sort by atom name (alphabetically)
#
# Run it something like:
#   sortPDB.pl <INFILE.pdb> > OUTFILE.pdb
# where INFILE.pdb and OUTFILE.pdb can be the same file.

my %atomIds = (
  " N  " => 0,
  " CA " => 1,
  " C  " => 2,
  " O  " => 3,
  " CB " => 4,
  " CG " => 5,
  " CD " => 6);

my $pdbin = shift;
my $pdbout = shift;

open(my $INF, "<", $pdbin) or die "Error: could not open input file: $!\n";
my $output = "";

my $prev_res = "";
my @lines = ();
my $atomnum = 1; # also fix the atom numbering along the way.
my $next_stuff = "";

while(<$INF>) {
  if (!/^ATOM|HETATM/) {
    if (@lines > 0) {
      $next_stuff .= $_;
    } else {
      $next_stuff = "";
      $output .= $_;
    }
    next;
  }

	my $line = $_;
	my $resi = substr($line, 22, 4);

  # If we've found a new residue.
  if ($prev_res != $resi) {
    if (@lines > 0) {
      # Sort the lines and add it to the output.
      $output .= outputSorted($atomnum, @lines);
      # Might need to output something *after* this residue.
      $output .= $next_stuff;
      $next_stuff = "";
    }
    # Increment count by the number of lines.
    $atomnum += @lines;
    # Reset lines
    @lines = ();
    $prev_res = $resi;
  }
  push @lines, $line;
}
# Check and see if there's something at the end.
if (@lines > 0) {
  $output .= outputSorted($atomnum, @lines);
}
$output .= $next_stuff;
close($INF);

open(my $OUTF, ">", $pdbout) or die "Error: could not open outfile! $!\n";
print $OUTF $output;
close($OUTF);

sub outputSorted {
  my $count = shift;
  my @lines = @_;
  my @sorted = sort byAtom @lines;
  my $str = "";
  for (@sorted) {
    my $line = $_;
    #print "line is [$line]\n";
    substr( $line, 5, 6 ) = sprintf("% 6d", $count++);
    $str .= $line;
  }
  return $str;
}
sub byAtom {
  my $atom_a = substr($a, 12, 4);
  my $atom_b = substr($b, 12, 4);

  my $sort_a = sortByAtomId($atom_a);
  my $sort_b = sortByAtomId($atom_b);
  if ($sort_a == $sort_b) {
    return $atom_a cmp $atom_b;
  }
  return $sort_a <=> $sort_b;
}
sub sortByAtomId {
  my $atomname = shift;
  if (defined $atomIds{$atomname}) {
    return $atomIds{$atomname};
  }
  # Otherwise, just return a number greater than any of these.
  return scalar keys %atomIds;
}

# Read the PDB in
