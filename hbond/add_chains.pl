#! /usr/bin/perl
################################
# Charm deletes the chain information. We'd like it back. This script will
# do that for us.
################################

my $pdb = shift;
my $chain = shift; # Add this chain
my $minres = shift; # ...starting at this residue number
my $maxres = shift; # ...and continuing until this residue number
# Will print to STDOUT

open(my $PDBIN, "<", $pdb) or die "Error: Could not open pdb input file [$pdb]: $!\n";
while(<$PDBIN>) {
  if (/^ATOM/) {
    my $resi = substr($_, 22, 4);

    # Alter the chain if needed.
    if ($resi >= $minres && $resi <= $maxres) {
      substr($_, 21, 1) = $chain;
    }
  }
  print $_;
}
